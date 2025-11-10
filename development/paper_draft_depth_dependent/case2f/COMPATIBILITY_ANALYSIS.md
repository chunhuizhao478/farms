# Compatibility Analysis: InitialStressStrainTPV26 & EffectiveBodyForceTPV26

## Executive Summary

✅ **Vertical Equilibrium**: PERFECT (residual = 0.00e+00 Pa/m)
⚠️ **Horizontal Stress (σ_xx)**: Small deviation observed
✅ **Other Components (σ_yy, σ_zz, σ_xy)**: Good agreement

## Equilibrium Verification Results

For `lambda_pp = 0.9`:

### Region 1 (0-6 km): Hydrostatic
- Effective body force: **16.37 kPa/m³** (constant)
- dσ'_zz/dz: **-16.37 kPa/m³**
- **Equilibrium residual: 0.00e+00 Pa/m** ✓

### Region 2 (6-8 km): Quadratic Transition
- Effective body force range: -102.65 to 25.78 kPa/m³
- dσ'_zz/dz range: -25.78 to 102.65 kPa/m³
- **Equilibrium residual: 0.00e+00 Pa/m** ✓

### Region 3 (>8 km): Scaled Lithostatic
- Effective body force: **2.62 kPa/m³** (10% effective stress)
- dσ'_zz/dz: **-2.62 kPa/m³**
- **Equilibrium residual: 0.00e+00 Pa/m** ✓

## Implementation Compatibility

### 1. Pore Pressure Calculation

Both classes use **identical** pore pressure formulations:

**Region 1 (0 to A)**: Hydrostatic
```
Pf = ρ_fluid · g · z
dPf/dz = ρ_fluid · g
```

**Region 2 (A to B)**: Quadratic transition
```
s = (z - A) / (B - A)
Pf = Pf_A + (Pf_B_target - Pf_A) · s²
dPf/dz = 2 · (Pf_B_target - Pf_A) / (B - A) · s
```

**Region 3 (> B)**: Scaled lithostatic
```
Pf = λ · ρ · g · z
dPf/dz = λ · ρ · g
```

✅ **Verified**: Implementations match exactly in both C++ classes.

### 2. Vertical Equilibrium

The vertical equilibrium equation:
```
dσ'_zz/dz + f_z = 0
```

Where:
- σ'_zz = σ_zz + Pf = -ρ·g·z + Pf
- dσ'_zz/dz = -ρ·g + dPf/dz
- f_z = ρ·g - dPf/dz (from EffectiveBodyForceTPV26)

Substituting:
```
(-ρ·g + dPf/dz) + (ρ·g - dPf/dz) = 0 ✓
```

✅ **Verified**: Perfect equilibrium achieved mathematically and numerically.

### 3. Horizontal Stress Calculation

Both classes use identical formulations:

**Total stress:**
```
σ_xx^total = Ω · (b_xx · (σ_zz + Pf) - Pf) + (1 - Ω) · σ_zz
σ_yy^total = Ω · (b_yy · (σ_zz + Pf) - Pf) + (1 - Ω) · σ_zz
σ_xy = Ω · (b_xy · (σ_zz + Pf))
```

**Effective stress:**
```
σ'_xx = σ_xx^total + Pf = [Ω · b_xx + (1 - Ω)] · (σ_zz + Pf)
σ'_yy = σ_yy^total + Pf = [Ω · b_yy + (1 - Ω)] · (σ_zz + Pf)
σ'_xy = σ_xy  (no change for shear)
```

✅ **Verified**: Implementations match exactly.

## Source of σ_xx Deviation

The small deviation in σ_xx (but not σ_yy, σ_zz, σ_xy) likely comes from one of these sources:

### 1. Eigenstrain Computation (Most Likely Cause)

The `ComputeDamageBreakageEigenstrainFromInitialStress` material:
1. Reads initial stress σ⁰ from functions
2. Solves a **nonlinear equation** for strain invariant ratio ξ using Newton iteration
3. Computes compliance tensor with damage-dependent parameters
4. Computes eigenstrain ε^e = -ε

**Key Issue**: The Newton iteration may not converge to exactly the same ξ value used in the Python script, especially if:
- Initial guess differs
- Convergence tolerance differs
- Damage parameters affect the solution path

Since all stress components depend on ξ, why only σ_xx deviates?
- The deviation might be amplified by the b_xx coefficient (0.926793)
- Boundary condition discretization on left/right faces may magnify errors
- The static solver may adjust σ_xx more than other components to satisfy global equilibrium

### 2. Boundary Condition Discretization

Stress boundary conditions are applied on left/right (x), front/back (y) faces:
- Left face: applies -σ_xx (normal = -x)
- Right face: applies +σ_xx (normal = +x)
- Front face: applies +σ_yy (normal = +y)
- Back face: applies -σ_yy (normal = -y)

**Possible issue**: If mesh nodes don't lie exactly at function evaluation points, interpolation errors could cause slight mismatches. This would affect boundary stresses more than interior stresses.

### 3. Static Solver Adjustment

The static solver finds displacement u such that:
```
K · u = F_ext - F_int
```

Where:
- K = stiffness matrix
- F_ext = external forces (boundary conditions + body force)
- F_int = internal forces (from eigenstrain and elastic response)

If the initial stress field is not **exactly** in discrete equilibrium (due to numerical errors), the solver will adjust it. The adjustment might be larger for σ_xx due to:
- Stronger coupling through b_xx coefficient
- Boundary effects on left/right faces
- Asymmetry in mesh or loading

### 4. Coordinate System and abs(z_coord)

Both classes use `abs(z_coord)` to get depth. If there are any nodes with z ≈ 0, numerical roundoff could cause inconsistencies. However, this would affect all stress components equally.

## Recommended Checks

### 1. Quantify the Deviation
Check the magnitude of σ_xx deviation:
```python
# In ParaView or Python postprocessing
deviation = (computed_sigma_xx - initial_sigma_xx) / initial_sigma_xx * 100
```

**Acceptable**: < 1% deviation
**Warning**: 1-5% deviation (check eigenstrain convergence)
**Error**: > 5% deviation (indicates implementation bug)

### 2. Check Spatial Distribution
Examine where the deviation occurs:
- If deviation is largest near boundaries → boundary condition issue
- If deviation is largest in transition zones (6-8 km, 15-20 km) → gradient discretization issue
- If deviation is uniform → eigenstrain computation issue

### 3. Verify Other Invariants
Check if the strain invariant ratio ξ matches:
```python
# From static solve output
xi_computed = static_solve_output['strain_invariant_ratio']

# From Python script
xi_analytical = I1 / sqrt(I2)

# Compare
max_xi_deviation = max(abs(xi_computed - xi_analytical))
```

### 4. Check Material Properties
Verify that damage-dependent parameters in eigenstrain computation match:
```
initial_damage = 0  (from input file line 229)
mu = shear_modulus_o + xi_o * damage * gamma_r
   = 3.204e10 + (-0.8) * 0 * gamma_r
   = 3.204e10  ✓ (matches input file)
```

## Expected Behavior

For a well-implemented static solve with overpressure:

1. **Vertical stress σ_zz**: Should match initial stress exactly or within < 0.1%
2. **Horizontal stresses σ_xx, σ_yy**: May have small deviations (< 1-2%) due to eigenstrain computation
3. **Shear stress σ_xy**: Should match initial stress very well (< 0.5%)
4. **Vertical displacement**: Should be small but non-zero (settling due to body force)
5. **Horizontal displacement**: Should be very small near boundaries where BCs are applied

## Fixes if Deviation is Too Large

### Fix 1: Use Simpler Eigenstrain Computation (if available)
If MOOSE has a `ComputeEigenstrainFromInitialStress` that uses simple linear elasticity instead of damage-dependent compliance, use that instead.

### Fix 2: Tighten Newton Solver Tolerance
In `ComputeDamageBreakageEigenstrainFromInitialStress.C`, the Newton loop uses:
```cpp
if (std::abs(F) < 1e-12)
  break; // converged
```

Try tightening this to `1e-14` or `1e-15`.

### Fix 3: Improve Initial Guess for ξ
Line 144 in eigenstrain material:
```cpp
Real xi = (std::abs(S1) > 1e-16 && std::abs(S2) > 1e-16) ? (S1 / std::sqrt(S2)) : _xi_o;
```

The initial guess could be improved by using the value from the Python script.

### Fix 4: Check λ and μ Consistency
Ensure that:
```
lambda_o (input file) = lambda_o (Python script) = 3.204e10
shear_modulus_o (input file) = shear_modulus_o (Python script) = 3.204e10
```

## Conclusion

The implementations of `InitialStressStrainTPV26` and `EffectiveBodyForceTPV26` are **mathematically compatible** and produce **perfect vertical equilibrium**. The small deviation in σ_xx is likely due to **numerical artifacts in the eigenstrain computation** using a Newton iteration for ξ.

As long as the deviation is < 1-2%, this is acceptable and will not significantly affect the dynamic simulation results. The deviation will be further smoothed out during the dynamic solve as the fault slips and damage evolves.

**Recommendation**: Document the deviation magnitude and proceed with the simulation if it's < 2%.
