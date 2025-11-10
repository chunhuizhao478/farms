# EffectiveBodyForceTPV26 Kernel Usage Guide

## Purpose

The `EffectiveBodyForceTPV26` kernel applies the effective body force that accounts for pore pressure gradients in overpressure models. It replaces the standard `BodyForce` kernel when using depth-dependent overpressure.

## Physics Background

The momentum equation with pore pressure is:
```
ρ ∂²u/∂t² = ∇·σ' + ρ_rock·g - ∇Pf
```

Where:
- σ' is effective stress (computed by MOOSE stress materials)
- ρ_rock·g is total body force
- ∇Pf is pore pressure gradient

The effective body force is:
```
f_eff = ρ_rock · g - ∇Pf
```

## Pore Pressure Gradients by Region

### Without Overpressure (Hydrostatic)
```
∇Pf = ρ_fluid · g = 9.8 kPa/m
f_eff = (ρ_rock - ρ_fluid) · g = 16.366 kPa/m³
```

### With Standard Overpressure (`use_overpressure = true`, `overpressure_loweffective = false`)

**Region 1 (0 to A)**: Hydrostatic
```
∇Pf = ρ_fluid · g
f_eff = (ρ_rock - ρ_fluid) · g
```

**Region 2 (A to B)**: Linear transition
```
∇Pf = g · (ρ_fluid + Δρ · (z - A) / (B - A))
f_eff = ρ_rock · g - ∇Pf  [spatially varying]
```

**Region 3 (> B)**: Fully lithostatic
```
∇Pf = ρ_rock · g
f_eff = 0  [zero effective stress]
```

### With Low Effective Stress Overpressure (`overpressure_loweffective = true`)

**Region 1 (0 to A)**: Hydrostatic (same as above)

**Region 2 (A to B)**: Quadratic transition
```
s = (z - A) / (B - A)
Pf_A = ρ_fluid · g · A
Pf_B_target = λ · ρ_rock · g · B
∇Pf = 2 · (Pf_B_target - Pf_A) / (B - A) · s
f_eff = ρ_rock · g - ∇Pf  [spatially varying]
```

**Region 3 (> B)**: Scaled lithostatic
```
∇Pf = λ · ρ_rock · g
f_eff = (1 - λ) · ρ_rock · g  [small positive effective stress]

For λ = 0.98:
f_eff = 0.02 · 2670 · 9.8 = 523 Pa/m³ (only 2% of overburden!)
```

## Input File Usage

### Example 1: Hydrostatic (Standard Case)

Replace the old body force:
```ini
# OLD - DON'T USE WITH OVERPRESSURE
[Kernels]
  [gravity_z]
    type = BodyForce
    variable = disp_z
    value = ${fparse -1 * density * gravity}
  []
[]
```

With the new kernel:
```ini
# NEW - Hydrostatic effective body force
[Kernels]
  [gravity_z]
    type = EffectiveBodyForceTPV26
    variable = disp_z
    fluid_density = 1000
    rock_density = 2670
    gravity = 9.8
    use_overpressure = false
  []
[]
```

### Example 2: Standard Overpressure (Full Lithostatic Below 8 km)

```ini
[Kernels]
  [gravity_z]
    type = EffectiveBodyForceTPV26
    variable = disp_z
    fluid_density = 1000
    rock_density = 2670
    gravity = 9.8
    use_overpressure = true
    overpressure_depth_A = 6000
    overpressure_depth_B = 8000
  []
[]
```

### Example 3: Low Effective Stress Overpressure (2% Effective Stress Below 8 km)

```ini
[Kernels]
  [gravity_z]
    type = EffectiveBodyForceTPV26
    variable = disp_z
    fluid_density = 1000
    rock_density = 2670
    gravity = 9.8
    use_overpressure = true
    overpressure_depth_A = 6000
    overpressure_depth_B = 8000
    overpressure_loweffective = true
    lambda_pp = 0.98    # 98% pore pressure, 2% effective stress
  []
[]
```

### Example 4: Complete Input File Section

```ini
[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
  volumetric_locking_correction = false
  order = FIRST
  family = LAGRANGE
[]

[Kernels]
  [SolidMechanics]
    # Automatically generates momentum equation kernels for x, y, z
  [../]

  [inertia_x]
    type = InertialForce
    variable = disp_x
  []
  [inertia_y]
    type = InertialForce
    variable = disp_y
  []
  [inertia_z]
    type = InertialForce
    variable = disp_z
  []

  # Effective body force accounting for pore pressure gradient
  [gravity_z]
    type = EffectiveBodyForceTPV26
    variable = disp_z
    fluid_density = 1000
    rock_density = 2670
    gravity = 9.8
    use_overpressure = true
    overpressure_depth_A = 6000
    overpressure_depth_B = 8000
    overpressure_loweffective = true
    lambda_pp = 0.98
  []
[]
```

## Parameter Summary

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `fluid_density` | Real | Required | Fluid density (kg/m³), typically 1000 |
| `rock_density` | Real | Required | Rock density (kg/m³), typically 2670 |
| `gravity` | Real | Required | Gravitational acceleration (m/s²), typically 9.8 |
| `use_overpressure` | bool | false | Enable overpressure model |
| `overpressure_depth_A` | Real | -1 | Transition start depth (m) |
| `overpressure_depth_B` | Real | -1 | Transition end depth (m) |
| `overpressure_loweffective` | bool | false | Enable low effective stress model |
| `lambda_pp` | Real | 0.9 | Pore pressure ratio (0 < λ ≤ 1) |

## Consistency Check

**IMPORTANT**: The parameters in `EffectiveBodyForceTPV26` kernel must match those in `InitialStressStrainTPV26` functions!

```ini
[Functions]
  # Initial stress uses these parameters
  [./ini_stress_xx]
    type = InitialStressStrainTPV26
    # ...
    overpressure_depth_A = 6000
    overpressure_depth_B = 8000
    overpressure_loweffective = true
    lambda_pp = 0.98
  [../]
[]

[Kernels]
  # Body force kernel MUST use same parameters
  [gravity_z]
    type = EffectiveBodyForceTPV26
    # ...
    overpressure_depth_A = 6000     # MUST MATCH
    overpressure_depth_B = 8000     # MUST MATCH
    overpressure_loweffective = true # MUST MATCH
    lambda_pp = 0.98                # MUST MATCH
  []
[]
```

## Physical Interpretation

### At 10 km depth with `lambda_pp = 0.98`:

**Total Vertical Stress:**
```
σ_zz = -ρ_rock · g · z = -2670 · 9.8 · 10000 = -262 MPa
```

**Pore Pressure:**
```
Pf = λ · ρ_rock · g · z = 0.98 · 2670 · 9.8 · 10000 = 256 MPa
```

**Effective Vertical Stress:**
```
σ'_zz = σ_zz + Pf = -262 + 256 = -5.2 MPa (only 2% of total!)
```

**Effective Body Force:**
```
f_eff = (1 - λ) · ρ_rock · g = 0.02 · 2670 · 9.8 = 523 Pa/m³
```

This creates an extremely weak effective confining pressure, facilitating damage and fault slip.

## Files Modified

- `/include/kernels/EffectiveBodyForceTPV26.h` (new)
- `/src/kernels/EffectiveBodyForceTPV26.C` (new)

## Related Classes

- `InitialStressStrainTPV26`: Computes initial stress/strain with same overpressure model
- Standard MOOSE `BodyForce` kernel: Use this for non-overpressure cases
