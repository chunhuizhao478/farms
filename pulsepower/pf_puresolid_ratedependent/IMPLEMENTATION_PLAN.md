# Implementation Plan: Rate-Dependent Phase Field Fracture

## Reference
Hofacker, M., & Miehe, C. (2013). A phase field model of dynamic fracture: Robust field updates for the analysis of complex crack patterns. *International Journal for Numerical Methods in Engineering*, 93(3), 276-301. DOI: 10.1002/nme.4387

---

## 1. Theoretical Background

### 1.1 Current Quasi-Static Formulation

The current phase field evolution equation (Eq. 40 in paper):

```
(G_c/l)[d - l^2 * Laplacian(d)] = 2(1-d) * H
```

where:
- `d` : phase field (damage variable, 0=intact, 1=broken)
- `G_c` : critical energy release rate (fracture toughness)
- `l` : regularization length scale
- `H` : history field containing maximum tensile strain energy

This is implemented via two kernels:
- **ADPFFDiffusion**: Gradient term `(grad(w), (2*G_c*l/c_0) * grad(d))`
- **ADPFFSource**: Source term `(w, d_psi/d_d)` from free energy derivative

### 1.2 Target Rate-Dependent Formulation

The viscous regularized equation (Eq. 46 in paper):

```
(G_c/l)[d - l^2 * Laplacian(d)] + eta * d_dot = 2(1-d) * H
```

**New term**: `eta * d_dot` (viscous crack resistance)

where:
- `eta` : viscosity parameter [Pa*s] or [N*s/m^2]
- `d_dot` : time derivative of phase field

### 1.3 Physical Interpretation

The viscous term provides:
1. **Time regularization**: Prevents instantaneous crack propagation
2. **Enhanced stability**: Improves numerical robustness for dynamic problems
3. **Rate-dependent behavior**: Crack propagation speed depends on loading rate
4. **Energy dissipation**: Additional viscous dissipation during fracture

---

## 2. Implementation Components

### 2.1 New Kernel: ADPFFViscousResistance

**Location**:
- Header: `include/kernels/phasefield/ADPFFViscousResistance.h`
- Source: `src/kernels/phasefield/ADPFFViscousResistance.C`

**Weak Form**:
```
integral_Omega (w * eta * d_dot) dV
  = integral_Omega (w * (eta/dt) * (d - d_old)) dV
```

Using backward Euler time discretization: `d_dot = (d - d_old) / dt`

**Residual Contribution**:
```cpp
return _eta[_qp] * (_u[_qp] - _d_old[_qp]) / _dt;
```

### 2.2 Input Parameters

| Parameter | Symbol | Units | Description |
|-----------|--------|-------|-------------|
| viscosity | eta | Pa*s (or N*s/mm^2) | Viscous resistance parameter |

---

## 3. Parameter Selection Guidelines

### 3.1 Reference Values from Paper (Section 5)

| Material | rho (kg/m^3) | E (GPa) | G_c (kN/m) | eta (Ns/mm^2) | dt (s) |
|----------|--------------|---------|------------|---------------|--------|
| Steel (Kalthoff) | 8000 | 190 | 22.17 | 1e-9 | ~1e-8 |
| PMMA (tension) | 1190 | 3.24 | 0.35 | 1e-10 | ~1e-9 |

### 3.2 Choosing Viscosity Parameter

**Dimensionless Analysis**:
The key dimensionless number is:
```
M = eta * v_c / G_c
```
where `v_c` is the characteristic crack velocity (typically Rayleigh wave speed).

**Practical Guidelines**:

1. **Rate-independent limit** (eta -> 0):
   - Set `eta = 0` or `eta ~ 1e-12 Ns/mm^2`
   - Recovers quasi-static phase field behavior

2. **Mild rate-dependence**:
   - `eta = 1e-10 to 1e-9 Ns/mm^2`
   - Slight stabilization, minimal effect on crack path

3. **Strong rate-dependence**:
   - `eta = 1e-8 to 1e-6 Ns/mm^2`
   - Significant rate effects, slower crack propagation
   - May suppress crack branching

### 3.3 Stability Criterion

For numerical stability, the viscous term should be comparable to the crack resistance:
```
eta / dt ~ G_c / l
```

This gives an estimate for the viscosity:
```
eta ~ (G_c / l) * dt
```

**Example for your case**:
- G_c = 100 J/m^2 = 100 N/m
- l = 2e-4 m
- dt = 5e-8 s
- eta ~ (100 / 2e-4) * 5e-8 = 2.5e-5 Pa*s = 2.5e-11 N*s/mm^2

Start with `eta = 1e-9 Ns/mm^2` and adjust based on results.

---

## 4. Algorithmic Implementation

### 4.1 Staggered Scheme (Table I in paper)

The operator-split algorithm in time step [t_n, t_{n+1}]:

1. **Update history field**:
   ```
   H = max(psi_0^+(strain(u_n)), H_n)
   ```

2. **Compute phase field** (minimization problem):
   ```
   d = argmin { integral [G_c*gamma(d,grad_d) + (eta/2*dt)*(d-d_n)^2 + (1-d)^2*H] dV }
   ```

3. **Compute displacement field** (elastodynamics):
   ```
   u = argmin { integral [psi(strain;d) + (rho*4/dt^2)*(0.5*|u|^2 - u_tilde_n*u) - gamma_bar*u] dV }
   ```

4. **Update history variables**: Store H, d, u, u_dot, u_ddot

### 4.2 Weak Form Implementation

**Total phase field residual**:
```
R = R_diffusion + R_source + R_viscous
```

where:
- `R_diffusion = (grad_w, 2*G_c*l/c_0 * grad_d)` [ADPFFDiffusion]
- `R_source = (w, d_psi/d_d)` [ADPFFSource]
- `R_viscous = (w, eta * (d - d_old) / dt)` [ADPFFViscousResistance] **NEW**

---

## 5. Input File Template

### 5.1 Material Definition

```ini
# Rate-dependent fracture parameters
l = 2e-4        # [m] regularization length
Gc_const = 100  # [J/m^2] fracture toughness
eta = 1e-9      # [Ns/mm^2] viscosity parameter

[Materials]
  [fracture_properties]
    type = ADGenericConstantMaterial
    prop_names = 'l Gc eta'
    prop_values = '${l} ${Gc_const} ${eta}'
  []
[]
```

### 5.2 Kernel Definition

```ini
[Kernels]
  [diff]
    type = ADPFFDiffusion
    variable = d
    fracture_toughness = Gc
    regularization_length = l
    normalization_constant = c0
  []
  [source]
    type = ADPFFSource
    variable = d
    free_energy = psi
  []
  [viscous]
    type = ADPFFViscousResistance
    variable = d
    viscosity = eta
  []
[]
```

---

## 6. Validation Strategy

### 6.1 Unit Tests

1. **Test 1**: Verify zero residual when d = d_old (no change)
2. **Test 2**: Verify linear scaling with viscosity parameter eta
3. **Test 3**: Verify correct time step scaling (1/dt dependence)
4. **Test 4**: Verify rate-independent limit (eta = 0)

### 6.2 Benchmark Problems

1. **Kalthoff-Winkler Test** (Section 5.1 of paper):
   - Impact velocities: 5, 16.5, 50 m/s
   - Expected: ~70 degree crack angle for brittle failure
   - Crack branching at higher velocities

2. **Single-edge notched tension** (Section 5.2 of paper):
   - Constant velocity loading
   - Compare crack paths for different loading rates

### 6.3 Convergence Studies

1. Mesh convergence: verify mesh-independent results
2. Time step convergence: verify dt-independent crack path
3. Viscosity sensitivity: quantify effect of eta on results

---

## 7. File Structure

```
farms_cdms/
├── include/kernels/phasefield/
│   └── ADPFFViscousResistance.h      # NEW
├── src/kernels/phasefield/
│   └── ADPFFViscousResistance.C      # NEW
├── pulsepower/
│   ├── pf_puresolid_ratedependent/
│   │   ├── IMPLEMENTATION_PLAN.md    # This document
│   │   ├── elasticity.i              # Elastodynamics input
│   │   └── fracture.i                # Phase field input (rate-dependent)
│   └── unit_tests/
│       └── kernels/
│           └── test_ADPFFViscousResistance.i  # Unit test
```

---

## 8. Expected Behavior

### 8.1 Effect of Viscosity on Crack Propagation

| eta (Ns/mm^2) | Behavior |
|---------------|----------|
| 0 | Rate-independent, instantaneous crack evolution |
| 1e-12 | Near rate-independent, minimal stabilization |
| 1e-9 | Mild rate-dependence, good stability |
| 1e-6 | Strong rate-dependence, slow crack propagation |

### 8.2 Comparison with Experiments

The Hofacker-Miehe model with rate-dependence should capture:
- Crack initiation delay under rapid loading
- Reduced crack branching with increasing viscosity
- Crack arrest phenomena
- Loading-rate dependent fracture toughness

---

## 9. References

1. Hofacker, M., & Miehe, C. (2013). A phase field model of dynamic fracture: Robust field updates for the analysis of complex crack patterns. IJNME, 93(3), 276-301.

2. Miehe, C., Hofacker, M., & Welschinger, F. (2010). A phase field model for rate-independent crack propagation: Robust algorithmic implementation based on operator splits. CMAME, 199(45-48), 2765-2778.

3. Miehe, C., Welschinger, F., & Hofacker, M. (2010). Thermodynamically consistent phase-field models of fracture: Variational principles and multi-field FE implementations. IJNME, 83(10), 1273-1311.

---

## 10. Implementation Checklist

- [ ] Create ADPFFViscousResistance.h header file
- [ ] Create ADPFFViscousResistance.C source file
- [ ] Register kernel in farmsApp
- [ ] Rebuild application (make -j8)
- [ ] Create unit test input file
- [ ] Run unit tests
- [ ] Create example rate-dependent input files
- [ ] Validate with Kalthoff-Winkler benchmark
- [ ] Document parameter sensitivity

---

*Document created: 2026-01-17*
*Author: Claude Code Assistant*
