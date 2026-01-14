# Large Deformation Implementation Notes

## Overview
This directory contains the Hencky-type finite deformation implementation for the three-point bending problem with phase field fracture.

## Files Created
1. **elasticity.i** - Main mechanics solver with large deformation framework
2. **fracture.i** - Phase field fracture solver (minimal changes from small deformation version)

## Key Modifications from Small Deformation Version

### elasticity.i Changes

#### 1. GlobalParams Section
- **Added**: `volumetric_locking_correction = true`
  - Prevents volumetric locking in nearly incompressible materials
  - Essential for accurate large deformation analysis

#### 2. Kernels Section
- **Modified**: All three `ADStressDivergenceTensors` kernels now have `use_displaced_mesh = true`
  - Enforces equilibrium on the deformed (current) configuration
  - Critical for geometric nonlinearity

#### 3. Materials Section - Block 1 (Phase Field Region)

**Replaced:**
- `ADComputeSmallStrain` → `ComputeDeformationGradient`
  - Computes deformation gradient F from displacement gradients

- `SmallDeformationIsotropicElasticity` → `HenckyIsotropicElasticity`
  - Uses logarithmic (Hencky) strain measure: ε = 0.5 * log(F^T * F)
  - Maintains same degradation and decomposition (spectral) as before
  - Outputs `psie_active` for phase field coupling

- `ComputeSmallDeformationStress` → `ComputeLargeDeformationStress`
  - Computes Cauchy stress from Mandel stress and deformation gradient

#### 4. Materials Section - Block 2 (Elastic Regions)

**Implemented:** Consistent large deformation framework without phase field
- `ComputeDeformationGradient` - Computes deformation gradient
- `NoDegradation` - Degradation function g = 1 (no damage)
- `HenckyIsotropicElasticity` - Same hyperelastic model, no decomposition needed
- `ComputeLargeDeformationStress` - Stress computation

### fracture.i Changes
- **No structural changes** - Phase field equations are independent of deformation framework
- `psie_active` is transferred from elasticity app and used identically

## Material Property Flow

### elasticity.i → fracture.i
1. `HenckyIsotropicElasticity` computes `psie_active` (line 189)
2. Transfer block sends `psie_active` to fracture app (lines 28-33)
3. Fracture app receives `psie_active` as AuxVariable (line 15)
4. Phase field energy `psi` uses `psie_active` (line 76)

### fracture.i → elasticity.i
1. Fracture app solves for damage field `d`
2. Transfer block sends `d` back to elasticity app (lines 22-27)
3. Degradation function `g(d)` modifies stress and energy

## Verification Checklist

### ✓ Consistency Verified
- [x] `psie_active` output from HenckyIsotropicElasticity
- [x] `psie_active` transfer to fracture app
- [x] `psie_active` used in phase field energy
- [x] Degradation functions identical in both apps
- [x] Spectral decomposition maintained
- [x] Volumetric locking correction enabled
- [x] Displaced mesh used in equilibrium

## Running the Simulation

### Command
```bash
cd /Users/chunhuizhao/projects/farms_cdms/threepointbending/code/case_wohole_largedeform
mpiexec -n <nprocs> ../../path/to/raccoon-opt -i elasticity.i
```

### Expected Behavior
1. Elasticity app initializes with large deformation framework
2. For each timestep:
   - Elasticity app solves for displacements using current damage field
   - Computes `psie_active` and transfers to fracture app
   - Fracture app updates damage field using `psie_active`
   - New damage field transferred back to elasticity app
3. Fixed-point iteration continues until convergence

## Comparison with Small Deformation

### To Validate Implementation
1. **Run with very small loads** (reduce loading by 100x)
   - Results should match small deformation version
   - Validates implementation correctness

2. **Run with actual loads** (current loading: -1e-4 * t)
   - Should show differences from small deformation
   - Captures geometric nonlinearity effects

### Expected Differences
- **Stress distribution**: Rotated with deformation
- **Load-displacement curve**: Nonlinear even before damage
- **Damage evolution**: May differ due to stress redistribution

## Solver Considerations

### Current Settings
- Solver: NEWTON with vinewtonrsls
- Preconditioner: hypre boomeramg
- Automatic scaling: enabled
- Fixed-point iterations: max 20
- Time step: 0.1

### If Convergence Issues Occur
1. **Reduce time step**: Change `dt` from 0.1 to 0.01 or 0.001
2. **Adjust nonlinear tolerances**: Relax `nl_rel_tol` to 1e-6
3. **Use direct solver**: Uncomment lu/superlu_dist options
4. **Add line search**: Remove `line_search = none` if added

## Output Variables

### Elasticity App (exodus file)
- `disp_x`, `disp_y`, `disp_z`: Displacement components
- `d`: Damage field (from fracture app)
- `stress`: Cauchy stress tensor
- `elastic_strain`: Logarithmic (Hencky) strain tensor
- `psie_active`: Active strain energy density
- `deformation_gradient`: F tensor (if added as AuxVariable)

### Fracture App
- `d`: Damage field (0=intact, 1=broken)
- `psie_active`: Transferred from elasticity app
- `psi`: Total phase field energy density

## Theory Notes

### Hencky Strain
The logarithmic (Hencky) strain is defined as:
```
E = 0.5 * log(F^T * F) = 0.5 * log(C)
```
where F is the deformation gradient and C is the right Cauchy-Green tensor.

### Spectral Decomposition
The elastic energy is split into positive (tensile) and negative (compressive) parts based on eigenvalues of the strain tensor. Only the positive part drives fracture:
```
psie_active = psie_positive
psie = g(d) * psie_positive + psie_negative
```

### Degradation Function
```
g(d) = (1-d)^2 * (1-eta) + eta
```
where `eta = 1e-6` prevents complete loss of stiffness.

## Known Limitations

1. **Large strains**: If strains exceed ~50%, may need Neo-Hookean or other hyperelastic models
2. **Contact**: Current implementation does not handle contact at crack faces
3. **Inertia**: Static analysis only (no dynamic effects)
4. **Temperature**: Isothermal analysis (no thermal effects)

## Mesh Requirements

The mesh file `../../meshfile/mesh_wohole_3d.msh` must:
- Be a valid GMSH format mesh
- Define boundaries: 2 (bottom_left_support), 3 (bottom_right_support), 4 (top_loading)
- Have suitable refinement for crack propagation (element size related to length scale `l`)

## Next Steps

1. **Validation run**: Test with reduced loads to verify against small deformation
2. **Full simulation**: Run with actual loads to observe large deformation effects
3. **Post-processing**: Extract force-displacement curves, damage patterns
4. **Parametric study**: Vary Gc, l, loading rate if needed

## References

- RACCOON tutorials: `raccoon/tutorials/large_deformation/`
- RACCOON J2 plasticity examples: `raccoon/tutorials/homogeneous_cube/Hencky_J2_*/`
- Hencky elasticity implementation: `raccoon/src/materials/large_deformation_models/HenckyIsotropicElasticity.C`
