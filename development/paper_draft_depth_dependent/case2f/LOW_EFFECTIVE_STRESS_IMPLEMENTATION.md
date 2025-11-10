# Low Effective Stress Overpressure Implementation

## Summary

Implemented the overpressure distribution from `plotsts_overpressure_loweffective.py` into `InitialStressStrainTPV26` C++ class.

## New Parameters

### 1. `overpressure_loweffective` (bool)
- **Default**: `false`
- **Description**: Flag to use low effective stress overpressure model (quadratic transition, lambda_pp scaling)
- **Requirement**: Must set `use_overpressure = true` when using this flag

### 2. `lambda_pp` (Real)
- **Default**: `0.9` (10% effective stress)
- **Valid Range**: (0, 1]
- **Description**: Pore pressure ratio controlling effective stress level below depth B
- **Physical Meaning**:
  - `lambda_pp = 0.90` → 10% effective stress (90% supported by fluid)
  - `lambda_pp = 0.95` → 5% effective stress (95% supported by fluid)
  - `lambda_pp = 0.98` → 2% effective stress (98% supported by fluid)

## Pore Pressure Model

### Three Regions:

1. **Region 1 (0 to overpressure_depth_A): Hydrostatic**
   ```
   Pf = ρ_fluid · g · z
   Gradient: 9.8 kPa/m (for ρ_fluid = 1000 kg/m³)
   ```

2. **Region 2 (overpressure_depth_A to overpressure_depth_B): Quadratic Transition**
   ```
   Pf_A = ρ_fluid · g · A
   Pf_B_target = λ · ρ_rock · g · B
   s = (z - A) / (B - A)
   Pf = Pf_A + (Pf_B_target - Pf_A) · s²
   ```
   - Gradual at start, steeper toward end
   - Smoothly connects hydrostatic to scaled lithostatic

3. **Region 3 (> overpressure_depth_B): Scaled Lithostatic**
   ```
   Pf = λ · ρ_rock · g · z
   Gradient: λ · 26.2 kPa/m (for ρ_rock = 2670 kg/m³)
   ```

## Effective Stress

Below depth B:
```
σ'_zz = σ_zz + Pf
      = -ρ_rock · g · z + λ · ρ_rock · g · z
      = -(1 - λ) · ρ_rock · g · z
```

**Effective stress is (1 - λ) fraction of total overburden**

## Example Usage in MOOSE Input File

```ini
[Functions]
  [./ini_stress_xx]
    type = InitialStressStrainTPV26
    i = 1
    j = 1
    lambda_o = 32.04e9
    shear_modulus_o = 32.04e9
    fluid_density = 1000
    rock_density = 2670
    gravity = 9.8
    bxx = 0.926793
    byy = 1.073206
    bxy = -0.8
    get_initial_stress = true
    use_tapering = true
    tapering_depth_A = 15000
    tapering_depth_B = 20000
    use_overpressure = true
    overpressure_depth_A = 6000
    overpressure_depth_B = 8000
    overpressure_loweffective = true    # Enable low effective stress
    lambda_pp = 0.98                     # 2% effective stress (98% pore pressure)
  [../]
[]
```

## Comparison with Standard Overpressure

### Standard Overpressure (`overpressure_loweffective = false`)
- **Region 2**: Linear transition with gradual density change
- **Region 3**: Fully lithostatic (λ = 1.0, zero effective stress)
- **Use Case**: Modeling fully lithostatic conditions

### Low Effective Stress (`overpressure_loweffective = true`)
- **Region 2**: Quadratic transition (more gradual)
- **Region 3**: Scaled lithostatic (λ < 1.0, positive effective stress)
- **Use Case**: Modeling weak but non-zero effective stress conditions

## Physical Interpretation

At 10 km depth with `lambda_pp = 0.98`:
- Total vertical stress: σ_zz = -262 MPa
- Pore pressure: Pf = 256 MPa
- Effective vertical stress: σ'_zz = -5.2 MPa (2% of total)

This creates an extremely overpressured environment where the rock skeleton only supports a small fraction of the overburden, facilitating fault slip and damage initiation.

## Modified Files

1. **Header**: `/include/functions/slipweakeningczm/InitialStressStrainTPV26.h`
   - Added `bool _overpressure_loweffective`
   - Added `Real _lambda_pp`

2. **Source**: `/src/functions/slipweakeningczm/InitialStressStrainTPV26.C`
   - Added parameter definitions
   - Added validation checks
   - Implemented quadratic transition and lambda_pp scaling in pore pressure calculation

## Notes

- The implementation exactly matches `plotsts_overpressure_loweffective.py`
- Quadratic transition provides smoother gradient change than linear
- Lambda_pp parameter allows flexible control of effective stress magnitude
- Works with all existing features (tapering, strain calculation, etc.)
