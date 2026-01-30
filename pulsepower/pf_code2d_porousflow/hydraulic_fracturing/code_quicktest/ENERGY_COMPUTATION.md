# Energy Computation in Hydraulic Fracturing Simulation

This document describes the energy balance computation implemented in the poroelastic hydraulic fracturing simulation.

## Overview

The simulation tracks energy balance using the first law of thermodynamics:

```
Input Energy = Stored Energy + Dissipated Energy
```

Or equivalently:

```
full_input_energy = full_energy
```

where `full_energy` contains both stored and dissipated components.

---

## Energy Balance Equation

### Total Input Energy (`full_input_energy`)

```
full_input_energy = -external_work - confinement_work + static_baseline - damping_work - fluid_injection_work
```

**Components:**

| Term | Description | Sign Convention |
|------|-------------|-----------------|
| `external_work` | Mechanical work by pulse loading on borehole | Negative = work done ON system |
| `confinement_work` | Work by confining stresses (top/bottom/left/right) | Negative = work done ON system |
| `damping_work` | Work absorbed by non-reflecting boundary dampers | Positive = energy leaving system |
| `fluid_injection_work` | Work by fluid injection at pore pressure BC | See detailed explanation below |
| `static_baseline` | Initial energy from static solve | Constant offset |

### Total Energy (`full_energy`)

```
full_energy = solid_kinetic + solid_elastic + solid_dissipated
            + fluid_kinetic + fluid_elastic + fluid_dissipated
```

---

## Solid Energy Terms

### 1. Solid Kinetic Energy

**Formula:**
```
E_kinetic = ∫ (1/2) ρ |v|² dV
```

**Implementation:**
- AuxVariable: `solid_kinetic_energy`
- AuxKernel: `KineticEnergyAux` using Newmark velocities (`vel_x`, `vel_y`, `vel_z`)
- Postprocessor: `solid_kinetic_energy_total` (ElementIntegralVariablePostprocessor)

### 2. Solid Elastic Energy

**Formula:**
```
E_elastic = ∫ ψ_e dV
```

where `ψ_e` is the strain energy density computed by the phase-field elasticity material.

**Implementation:**
- Material property: `psie` (from `NDSmallDeformationIsotropicElasticity`)
- Postprocessor: `solid_elastic_energy_dynamic` (ElementIntegralMaterialProperty)
- Note: Uses spectral decomposition for tension/compression split

### 3. Solid Dissipated Energy (Fracture)

**Formula (AT1 model):**
```
E_fracture = ∫ (Gc/c0) * (α/l + l|∇d|²) dV
```

where:
- `Gc` = critical energy release rate
- `c0` = normalization constant (8/3 for AT1)
- `l` = regularization length
- `α` = crack geometric function
- `d` = phase field (damage)

**Implementation:**
- Material: `CrackDissipatedEnergyDensity` in fracture subapp
- Postprocessor: `dissipated_energy_dynamic` (transferred from subapp)
- Total: `solid_dissipated_energy_total = dissipated_energy_dynamic + static_baseline`

---

## Fluid Energy Terms

### 4. Fluid Kinetic Energy

**Formula:**
```
E_fluid_kinetic = ∫ (1/2) ρ_f |q|² dV
```

where `q` is the Darcy velocity.

**Implementation:**
- AuxKernel: ParsedAux computing `0.5 * ρ_f * (qx² + qy² + qz²)`
- Postprocessor: `fluid_kinetic_energy_total`

**Note:** This is typically small compared to other terms.

### 5. Fluid Compression Energy

**Formula:**
```
E_fluid_compression = ∫ (1/2M) p² dV
```

where:
- `M` = Biot modulus (damage-dependent: `PorousFlow_constant_biot_modulus_qp`)
- `p` = pore pressure

This is the thermodynamically correct stored energy for fluid compression in porous media.
It represents the energy stored due to compressing the pore fluid.

**Implementation:**
- AuxVariable: `fluid_compression_energy`
- AuxKernel: ParsedAux computing `0.5 / M * p²`
- Postprocessor: `fluid_compression_energy_total`

### 6. Fluid Dissipated Energy (Darcy Viscous Dissipation)

**Formula:**
```
D_viscous = ∫∫ (μ/k) |q|² dV dt = ∫∫ q · ∇p dV dt
```

This is the energy dissipated by viscous friction during fluid flow through porous media.

**Implementation:**
- AuxKernel: ParsedAux computing `(μ/k) * |q|²`
- Uses local permeability from `effective_perm00_aux` (unified across all blocks)
- Postprocessor: `darcy_viscous_power` (instantaneous rate, units: W)
- Time-integrated: `darcy_viscous_dissipation_total` using `TimeIntegratedPostprocessor` (units: J)

**IMPORTANT:** Must use `TimeIntegratedPostprocessor` (not `CumulativeValuePostprocessor`) because
`darcy_viscous_power` is a rate [W = J/s]. Time integration computes ∫ power × dt to get energy [J].

**Total Fluid Dissipation:**
```
fluid_dissipated_energy_total = darcy_viscous_dissipation_total
```

**Note:** The previous "fluid elastic dissipation" term (property evolution effects) was removed because:
1. The term `0.5*α*(-p)*ε_v` has no clear thermodynamic basis in standard Biot poroelasticity
2. The coupling term `-αpε_v` appears in the stress relation, not as a separate stored energy
3. Property evolution effects from damage are already captured in the fracture dissipation energy

---

## Input Energy Terms

### External Work (Pulse Loading)

**Boundaries:** `hole`, `hole_fracture`

**Formula:**
```
W_external = ∫∫ f · u̇ dA dt
```

where `f` is the boundary force (from Pressure BC) and `u̇` is velocity.

**Implementation:**
- `FarmsExternalWork` postprocessor
- Forces saved via `save_in_disp_x/y` in Pressure BC
- Combined: `external_work = external_work_hole + external_work_hole_fractures`

### Confinement Work

**Boundaries:** `top`, `bottom`, `left`, `right`

Same formulation as external work but with constant confining pressure.

### Damping Work

**Boundaries:** `top`, `bottom`, `left`, `right`

Work absorbed by non-reflecting dashpot boundary conditions (`FarmsNonReflectDashpotBC`).

---

## Fluid Injection Work

### When to Include

| Scenario | Include in Input? | Reason |
|----------|-------------------|--------|
| No pore pressure BC | No | Internal redistribution only |
| Dirichlet pore pressure BC | **Yes** | External fluid source/sink |

### Computation

**Formula:**
```
W_injection = ∫∫ p (q · n) dA dt
```

where:
- `p` = pore pressure at boundary
- `q` = Darcy velocity
- `n` = outward normal (radially outward from origin for hole_fracture)

**Implementation:**

1. **Compute radial flux:**
   ```
   darcy_flux_normal = (qx * x + qy * y) / √(x² + y²)
   ```

2. **Compute power:**
   ```
   fluid_injection_power = ∫ p * darcy_flux_normal dA
   ```
   (SideIntegralVariablePostprocessor over `hole_fracture`, units: W)

3. **Time-integrate and negate:**
   ```
   fluid_injection_work = -TimeIntegrated(fluid_injection_power)
   ```
   **IMPORTANT:** Must use `TimeIntegratedPostprocessor` (not `CumulativeValuePostprocessor`) because
   `fluid_injection_power` is a rate [W = J/s]. Time integration computes ∫ power × dt to get energy [J].

**Sign Convention:**
- Positive `darcy_flux_normal` = flow INTO domain (injection)
- Positive `fluid_injection_work` = energy input to system

---

## Static Baseline Values

The simulation is initialized from a static solve. Baseline energy values are obtained from the static solve CSV output:

```
solid_dissipated_energy_total_static = <from static solve>
full_input_energy_static = <from static solve>
```

**Important:** After running the static solve, update these values in the dynamic input file header.

The static solve now computes:
- `fluid_compression_energy_total` using `(1/2M)p²` (thermodynamically correct)
- `solid_elastic_energy_total` from strain energy density
- `solid_dissipated_energy_total` from fracture energy

The dynamic solve uses these static baseline values to properly account for the initial energy state.

---

## Postprocessor Summary

### Input Energy Components

| Postprocessor | Description |
|---------------|-------------|
| `external_work` | Pulse loading work (hole + hole_fracture) |
| `confinement_work` | Confining pressure work (all sides) |
| `damping_work` | Dashpot absorbed energy |
| `fluid_injection_work` | Pore pressure BC work (if enabled) |
| `full_input_energy` | Total input energy |

### Stored/Dissipated Energy Components

| Postprocessor | Description |
|---------------|-------------|
| `solid_kinetic_energy_total` | Solid kinetic energy |
| `solid_elastic_energy_total` | Solid strain energy |
| `solid_dissipated_energy_total` | Fracture dissipation |
| `fluid_kinetic_energy_total` | Fluid kinetic energy |
| `fluid_compression_energy_total` | Fluid compression energy (1/2M)p² |
| `fluid_dissipated_energy_total` | Darcy viscous dissipation |
| `darcy_viscous_dissipation_total` | Darcy flow viscous dissipation (TimeIntegratedPostprocessor) |
| `full_energy` | Total energy (stored + dissipated) |

---

## CSV Output

The following quantities are output to CSV:

```
full_energy, full_input_energy,
solid_elastic_energy_total, solid_kinetic_energy_total, solid_dissipated_energy_total,
fluid_compression_energy_total, fluid_kinetic_energy_total, fluid_dissipated_energy_total,
darcy_viscous_dissipation_total, darcy_viscous_power,
damping_work, confinement_work, external_work,
fluid_injection_work, fluid_injection_power
```

---

## Energy Balance Verification

For a well-posed simulation, verify:

```
|full_input_energy - full_energy| / full_input_energy < tolerance
```

Common causes of imbalance:
1. Missing Darcy viscous dissipation (fixed in current implementation)
2. Incorrect sign conventions
3. Missing boundary work terms
4. Time integration errors
5. **Using `CumulativeValuePostprocessor` instead of `TimeIntegratedPostprocessor` for rate quantities**
   - `CumulativeValuePostprocessor` just sums: ∑(value_i)
   - `TimeIntegratedPostprocessor` integrates: ∫ value × dt = ∑(value_i × dt_i)
   - For power/rate quantities [W = J/s], must use `TimeIntegratedPostprocessor` to get energy [J]

---

## File Structure

- **Static solve:** `static_solve_quicktest/elasticityhf_static.i` - Computes initial equilibrium state
  - Outputs baseline energy values to CSV
  - Uses `fluid_compression_energy_total` = (1/2M)p²
- **Dynamic solve:** `code_quicktest/elasticityhf_wpulsep.i` - Main dynamic simulation
  - Contains all energy postprocessors
  - Uses static baseline values from static solve
- **Fracture subapp:** `fracturehf_*.i` - Computes fracture dissipation energy
- **Transfer:** Fracture energy transferred via `MultiAppPostprocessorTransfer`

### Workflow

1. Run static solve to get baseline energies
2. Update `full_input_energy_static` and `solid_dissipated_energy_total_static` in dynamic input
3. Run dynamic solve

---

## References

- Biot poroelasticity theory
- Phase-field fracture (AT1 model)
- Darcy flow energy dissipation
- CMAME paper (Appendix A for enhanced history energy)
