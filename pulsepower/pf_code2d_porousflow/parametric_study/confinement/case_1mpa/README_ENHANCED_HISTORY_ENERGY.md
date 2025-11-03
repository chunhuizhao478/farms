# Enhanced History Energy Implementation for Phase Field Damage

## Overview

This implementation updates the phase field damage model to use the enhanced maximum history energy formulation described in **Appendix A** of the CMAME paper "Numerical Simulation of Pulse Power-Induced Fracturing in Saturated Rock Sample".

## Mathematical Formulation

### Previous Implementation (Equation A.13)
```
H⁺ = max_t (2·ψₑ⁺)
```
Where:
- ψₑ⁺ is the positive (tensile) elastic energy density

### New Implementation (Equation A.12)
```
H⁺ = max_t ( ⟨2·ψₑ⁺ + p²·[(φₒ - 1)·(1/Kf - 1/Ks) - K/Ks²]⟩⁺ )
```

Where:
- ψₑ⁺ = positive elastic energy density
- p = pore pressure
- φₒ = initial porosity
- Kf = fluid bulk modulus
- Ks = grain bulk modulus
- K = solid bulk modulus
- ⟨·⟩⁺ = Macaulay bracket (positive part)

### Physical Meaning

The additional pressure-dependent term accounts for:
1. **Poroelastic property evolution** - How damage affects biot coefficient, porosity, and biot modulus
2. **Fluid-solid coupling** - The influence of pore pressure on fracture propagation
3. **Enhanced energy criterion** - More accurate representation of the energy available to drive damage

## Implementation Details

### New Source Files

#### 1. ElkPorousFlowHistoryEnergyEnhanced.h
- **Location**: `include/materials/porousflowmatprops/ElkPorousFlowHistoryEnergyEnhanced.h`
- **Purpose**: Header file defining the material class

#### 2. ElkPorousFlowHistoryEnergyEnhanced.C
- **Location**: `src/materials/porousflowmatprops/ElkPorousFlowHistoryEnergyEnhanced.C`
- **Purpose**: Implementation of the enhanced history energy calculation
- **Key Features**:
  - Computes pressure-dependent coefficient from poroelastic properties
  - Applies Macaulay bracket for thermodynamic consistency
  - Enforces irreversibility through max operation over time
  - Uses stateful material properties to track history

### Modified Input Files

#### elasticity.i

**Changes made:**

1. **Added AuxVariable** (line ~246-249):
```
[psie_active_enhanced]
  order = CONSTANT
  family = MONOMIAL
[]
```

2. **Added AuxKernel** (line ~367-372):
```
[psie_active_enhanced_aux]
  type = MaterialRealAux
  variable = psie_active_enhanced
  property = psie_active_enhanced
  execute_on = 'TIMESTEP_END'
[]
```

3. **Added Material Block** (line ~575-584):
```
[history_energy_enhanced]
  type = ElkPorousFlowHistoryEnergyEnhanced
  psie_active = psie_active
  pore_pressure = pp
  initial_porosity = ${porosity}
  fluid_bulk_modulus = ${fluid_bulk_modulus}
  grain_bulk_modulus = ${grain_bulk_modulus}
  bulk_modulus = K
  psie_active_enhanced = psie_active_enhanced
[]
```

4. **Modified Transfer** (line ~68):
Changed from:
```
source_variable = 'psie_active mesh_size'
```
To:
```
source_variable = 'psie_active_enhanced mesh_size'
```

#### fracture.i
**No changes needed** - The file already receives `psie_active` through the transfer and uses it correctly in the free energy expression (line 144).

## Usage

### Parameters Required

The new material requires the following parameters (all available in elasticity.i):

| Parameter | Variable in Input | Example Value |
|-----------|------------------|---------------|
| Initial porosity | `${porosity}` | 0.008 |
| Fluid bulk modulus | `${fluid_bulk_modulus}` | 1e9 Pa |
| Grain bulk modulus | `${grain_bulk_modulus}` | Derived from biot coefficient |
| Solid bulk modulus | `K` (material property) | ~31.4 GPa |
| Pore pressure | `pp` (variable) | Dynamic during simulation |
| Positive elastic energy | `psie_active` (material property) | Computed by elasticity material |

### Compilation

To compile the updated code:

```bash
cd /path/to/farms_cdms
make -j8
```

### Running the Simulation

```bash
cd pulsepower/pf_code2d_porousflow/parametric_study/confinement/case_1mpa
mpirun -n 8 ../../../../../farms-opt -i elasticity.i
```

## Verification and Testing

### Expected Behavior

1. **Energy Growth**: The enhanced history energy should be >= the basic version due to the positive pressure contribution
2. **Irreversibility**: H⁺ should be monotonically non-decreasing in time
3. **Pressure Dependence**: Higher pore pressures should lead to larger driving forces for damage
4. **Initial Conditions**: At t=0 with p=0, should match the basic formulation

### Debugging Tips

1. **Check Material Property**: Add output for `psie_active_enhanced` in the exodus file
2. **Verify Pressure Values**: Ensure pore pressure `pp` is properly coupled
3. **Compare with Basic**: Run both versions and compare damage evolution

### Common Issues

| Issue | Cause | Solution |
|-------|-------|----------|
| Compilation error | New files not registered | Check that `registerMooseObject` is present |
| Zero enhanced energy | Missing aux kernel | Verify AuxKernel extracts material property |
| Transfer failure | Variable name mismatch | Check transfer block variable names |
| Negative coefficient | Incorrect bulk moduli | Verify grain_bulk_modulus > K |

## Physical Interpretation

### Role of the Pressure Term

The pressure coefficient:
```
coeff = (φₒ - 1)·(1/Kf - 1/Ks) - K/Ks²
```

Typically:
- **(φₒ - 1)** is negative (porosity < 1)
- **(1/Kf - 1/Ks)** is positive (fluid more compressible)
- **K/Ks²** is positive but small

The net coefficient is usually **negative**, meaning:
- Positive pore pressures **increase** the driving force for damage
- This matches physical intuition: pressurized fluid promotes fracture

### Comparison with Previous Model

| Aspect | Previous (A.13) | Enhanced (A.12) |
|--------|----------------|-----------------|
| Energy source | Elastic strain only | Elastic + pore pressure |
| Coupling | Indirect (through stress) | Direct (in history variable) |
| Damage threshold | Purely mechanical | Hydro-mechanical |
| Physical realism | Approximate | More accurate |

## References

1. **CMAME Paper**: "Numerical Simulation of Pulse Power-Induced Fracturing in Saturated Rock Sample"
   - See **Appendix A**, equations (A.11) - (A.12)

2. **Related Work**:
   - Miehe et al. (2016) - Phase field fracture in hydro-poro-elasticity
   - Yu et al. (2024) - Hydraulic fracture in poroelastic media

## Future Extensions

Possible improvements to this implementation:

1. **Anisotropic permeability**: Account for directional permeability evolution
2. **Temperature coupling**: Include thermal effects on fluid properties
3. **Plastic effects**: Extend to elastoplastic-damage models
4. **Higher-order terms**: Include damage-rate dependent terms

## Contact

For questions about this implementation:
- **Author**: Chunhui Zhao
- **Institution**: University of Illinois Urbana-Champaign
- **Paper**: CMAME Pulsed Power Fracturing Study

---

**Last Updated**: 2025-11-02
**Code Version**: Enhanced history energy with poroelastic coupling
