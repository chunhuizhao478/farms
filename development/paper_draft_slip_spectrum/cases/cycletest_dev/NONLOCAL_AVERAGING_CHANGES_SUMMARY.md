# Nonlocal Averaging Modification Summary

## Overview
Modified the code to perform nonlocal averaging on **damage (α)** and **breakage (B)** fields instead of strain_invariant_ratio (ξ), and pass these nonlocal values back to the main app.

---

## New Source Files Created (Copies - Not In-Place Modifications)

### 1. Main App Material
- **Header**: `include/materials/cdbm_lagrangian/DiffusedDamageBreakageMaterialMainAppNonlocal.h`
- **Source**: `src/materials/cdbm_lagrangian/DiffusedDamageBreakageMaterialMainAppNonlocal.C`
- **Key Changes**:
  - Receives `alpha_damagedvar_nonlocal_aux` and `B_damagedvar_nonlocal_aux` from sub app
  - Uses nonlocal averaged values to compute material properties

### 2. Sub App Material
- **Header**: `include/materials/cdbm_lagrangian/DiffusedDamageBreakageMaterialSubAppNonlocal.h`
- **Source**: `src/materials/cdbm_lagrangian/DiffusedDamageBreakageMaterialSubAppNonlocal.C`
- **Key Changes**:
  - Couples local damage/breakage variables: `alpha_damagedvar_sub` and `B_damagedvar_sub`
  - Receives two `ElkRadialAverage` UserObjects: `average_alpha_UO` and `average_B_UO`
  - Declares nonlocal material properties: `alpha_nonlocal` and `B_nonlocal`
  - Performs radial averaging on damage and breakage fields

---

## Modified Input Files (In-Place)

### Main App: `unified_solve_main_eps1em7_cd10_cycle.i`

#### Changes:
1. **Materials** (line 374): `type = DiffusedDamageBreakageMaterialMainAppNonlocal`
2. **Removed**: Old nonlocal averaging on strain_invariant_ratio (lines 410-429)
3. **Transfers** (lines 755-770):
   - Pull: `alpha_nonlocal_sub`, `B_nonlocal_sub` from sub app
   - Push: `xi_aux` (not `nonlocal_xi`) to sub app
4. **Outputs** (line 553): Note that `alpha_damagedvar_aux` and `B_damagedvar_aux` now contain NONLOCAL values

### Sub App: `dynamic_solve_sub_eps1em7_cd10_cycle.i`

#### Changes:
1. **AuxVariables** (lines 132-148): Added output variables for comparison
2. **AuxKernels** (lines 240-266): Extract and copy variables for comparison
3. **Materials** (lines 271-303):
   - Changed to `DiffusedDamageBreakageMaterialSubAppNonlocal`
   - Added `ParsedMaterial` to convert variables to material properties
4. **UserObjects** (lines 326-343): Two `ElkRadialAverage` objects for α and B
5. **Outputs** (lines 392-402): Export comparison variables

---

## Output Variables for Comparison

### Sub App Outputs (for detailed comparison):

| Variable Name | Description | Purpose |
|---------------|-------------|---------|
| `alpha_damagedvar_sub` | Local damage (evolved) | **BEFORE averaging** |
| `B_damagedvar_sub` | Local breakage (evolved) | **BEFORE averaging** |
| `alpha_local_output` | Copy of local damage | **BEFORE averaging** (for easier identification) |
| `B_local_output` | Copy of local breakage | **BEFORE averaging** (for easier identification) |
| `alpha_nonlocal_sub` | Nonlocal averaged damage | **AFTER averaging** |
| `B_nonlocal_sub` | Nonlocal averaged breakage | **AFTER averaging** |

### Main App Outputs:

| Variable Name | Description | Note |
|---------------|-------------|------|
| `alpha_damagedvar_aux` | Damage field in main app | **Contains NONLOCAL averaged values** from sub app |
| `B_damagedvar_aux` | Breakage field in main app | **Contains NONLOCAL averaged values** from sub app |

---

## Data Flow Diagram

### Before (Old Implementation):
```
Main App:
  ├─ Compute strain_invariant_ratio (ξ)
  ├─ Nonlocal averaging: ξ → ξ_nonlocal
  └─ Send ξ_nonlocal to Sub App
        ↓
Sub App:
  ├─ Receive ξ_nonlocal
  ├─ Damage/breakage evolution using ξ_nonlocal
  └─ Send α, B (local) back to Main App
```

### After (New Implementation):
```
Main App:
  ├─ Compute strain_invariant_ratio (ξ)
  └─ Send ξ (local) to Sub App
        ↓
Sub App:
  ├─ Receive ξ (local)
  ├─ Damage/breakage evolution → α, B (local)
  ├─ Nonlocal averaging: α → α_nonlocal, B → B_nonlocal
  └─ Send α_nonlocal, B_nonlocal back to Main App
        ↓
Main App:
  └─ Use α_nonlocal, B_nonlocal in constitutive model
```

---

## Nonlocal Averaging Parameters (Consistent)

| Parameter | Value | Location |
|-----------|-------|----------|
| `length_scale` | 200 | Sub app UserObjects |
| `radius` | 400 | Sub app UserObjects |
| `weights` | BAZANT | Sub app UserObjects |
| `execute_on` | TIMESTEP_END | Sub app UserObjects |

---

## Material Parameters (Verified Consistent)

| Parameter | Main App | Sub App |
|-----------|----------|---------|
| `lambda_o` | 32.04e9 | 32.04e9 |
| `shear_modulus_o` | 32.04e9 | 32.04e9 |
| `xi_0` | -0.8 | -0.8 |
| `xi_d` | -0.9 | -0.9 |
| `C_g` | 1e-10 | (strain rate dependent) |
| `m1` | 10 | - |
| `m2` | 1 | - |
| `chi` | 0.8 | - |

---

## How to Compare Results

### In Paraview or Visualization Tool:

1. **Open Sub App output** (`dynamic_solve_sub_eps1em7_cd10_cycle_sub_app0.e`)

2. **Compare damage before/after averaging**:
   - Plot `alpha_local_output` (or `alpha_damagedvar_sub`) - LOCAL
   - Plot `alpha_nonlocal_sub` - NONLOCAL
   - Difference: `alpha_nonlocal_sub - alpha_local_output`

3. **Compare breakage before/after averaging**:
   - Plot `B_local_output` (or `B_damagedvar_sub`) - LOCAL
   - Plot `B_nonlocal_sub` - NONLOCAL
   - Difference: `B_nonlocal_sub - B_local_output`

4. **Verify in Main App**:
   - Open main app output (`unified_solve_main_eps1em7_cd10_cycle_out.e`)
   - Plot `alpha_damagedvar_aux` - should match `alpha_nonlocal_sub` from sub app
   - Plot `B_damagedvar_aux` - should match `B_nonlocal_sub` from sub app

---

## Expected Physical Effects

Nonlocal averaging on damage/breakage will:
- **Smooth sharp damage gradients** → more diffuse damage zones
- **Reduce mesh dependency** → more objective results
- **Regularize localization** → prevent zero-width shear bands
- **Preserve energy dissipation** → weighted averaging conserves integrated quantities

The smoothing effect should be visible by comparing:
- Sharp features in `alpha_local_output` vs. smooth `alpha_nonlocal_sub`
- Sharp features in `B_local_output` vs. smooth `B_nonlocal_sub`

---

## Compilation Instructions

After creating the new source files, recompile:

```bash
cd /Users/chunhuizhao/projects/farms_cdbm_implicit
make -j8
```

The new materials will be registered as:
- `DiffusedDamageBreakageMaterialMainAppNonlocal`
- `DiffusedDamageBreakageMaterialSubAppNonlocal`

---

## Testing Checklist

- [ ] Code compiles without errors
- [ ] Main app runs and converges
- [ ] Sub app runs and converges
- [ ] Transfer of nonlocal damage/breakage succeeds
- [ ] Output files contain all comparison variables
- [ ] `alpha_nonlocal_sub` is smoother than `alpha_local_output`
- [ ] `B_nonlocal_sub` is smoother than `B_local_output`
- [ ] Main app `alpha_damagedvar_aux` matches sub app `alpha_nonlocal_sub`
- [ ] Results are physically reasonable
- [ ] Nonlocal averaging reduces mesh sensitivity

---

## Contact / Questions

If you encounter issues:
1. Check that all new source files are in correct directories
2. Verify compilation completed successfully
3. Check input file syntax for typos
4. Verify variable names match in transfers
5. Check that ParsedMaterial expressions are correct

---

*Last updated: [Current Date]*
*Modification author: Claude Code*
