# Configurational Force Validation Tests

## Overview

These tests validate the `ConfigurationalForceUserObject` implementation, which computes configurational forces at nodes using the Eshelby stress tensor.

## Theory

For small deformation:
- **Eshelby stress**: Σ = ΨI - σ
- **Configurational force**: F_CNF = ∫ Σ · ∇N dV
- **Physical meaning**: Magnitude quantifies crack driving force, direction indicates propagation path

## Test Cases

### 1. `simple_tension_test.i`

**Purpose**: Validate basic functionality on a simple geometry

**Setup**:
- 3D rectangular bar (0.5 × 0.5 × 1.0 m)
- Material: Linear elastic (E = 10 GPa, ν = 0.3)
- Loading: Uniaxial tension in z-direction (0.5% strain)
- No actual crack - tests the Eshelby stress computation

**Expected Results**:
1. Energy density increases with loading
2. Configurational forces are non-zero
3. Force computation is stable and converges
4. Eshelby stress Σ = ΨI - σ is computed correctly

**Validation**:
```bash
cd test/tests/configurational_force
../../../farms-opt -i simple_tension_test.i
```

Check output:
- `force_mag` should be > 0
- `elastic_energy_aux` should increase with time
- No NaN or Inf values

### 2. `penny_crack_mode1.i` (Future)

**Purpose**: Validate against analytical solution for penny-shaped crack

**Setup**:
- Penny-shaped crack in 3D cube
- Mode I loading (pure tension)
- Compare with analytical K_I

**Expected Results**:
- Configurational force primarily in normal direction
- Magnitude scales with stress intensity factor
- Symmetry: F_x ≈ F_y ≈ 0 at crack center

## Analytical Validation

For a penny-shaped crack of radius `a` under remote tension σ_∞:

**Stress Intensity Factor**:
```
K_I = (2/π) * σ_∞ * √(πa)
```

**Energy Release Rate**:
```
G = K_I² / E'   where E' = E/(1-ν²) for plane strain
```

**Configurational Force**:
```
|F_CNF| ≈ G * (crack front length)
```

## Usage with Your Damage-Breakage Model

To use with `ComputeDamageBreakageStress3DSlipWeakeningNonlocal`:

```
[Materials]
  [stress]
    type = ComputeDamageBreakageStress3DSlipWeakeningNonlocal
    # ... your parameters ...
  []
[]

[UserObjects]
  [config_force]
    type = ConfigurationalForceUserObject
    energy_density = total_energy_density  # From line 324
    stress = sts_total                     # From line 310
    use_displacement_gradient = false
    execute_on = 'TIMESTEP_END'
  []
[]
```

## Visualization in Paraview

1. Load `output.e`
2. Apply **Glyph** filter:
   - Glyph Type: Arrow
   - Scalars: force_mag
   - Vectors: (force_x, force_y, force_z)
3. Color by `force_mag`

## References

- Santarossa et al. (2025) "Configurational forces explain echelon cracks in soft materials"
- Eshelby (1951) "The force on an elastic singularity"
- Rice (1968) "A path independent integral and the approximate analysis of strain concentration"
