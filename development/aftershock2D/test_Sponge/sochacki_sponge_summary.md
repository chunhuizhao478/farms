# Sochacki Sponge Boundary Summary

## Governing Idea

The Sochacki sponge extends the elastic wave equation with a spatially varying attenuation term:

```
d2u/dt2 + 2*A(x,y)*du/dt = v^2*(d2u/dx2 + d2u/dy2)
```

Here `A(x,y)` is zero in the physical domain and grows inside the sponge layer so outgoing waves are dissipated before reaching the exterior boundary. The factor of `2` follows the original Sochacki et al. (1987) formulation.

## Damping Profiles Implemented

The attenuation coefficient is built from a normalized distance `xi in [0,1]` between the interior edge of the sponge and the mesh boundary. Users can select one of the canonical ramps:

| Profile      | Expression (with S_max = sponge_s_max)                     |
|--------------|------------------------------------------------------------|
| linear       | `A = S_max * xi`                                           |
| exponent     | `A = S_max * xi^b` (`sponge_exponent_power = b`)           |
| cubic        | `A = S_max * xi^3`                                         |
| exponential  | `A = S_max * (exp(k*xi) - 1)/(exp(k) - 1)` (`k = sponge_exp_rate`) |
| gaussian     | `A = S_max * (1 - exp(-g*xi^2))/(1 - exp(-g))` (`g = sponge_gaussian_rate`) |

These match the dampers described in the "Comparison of artificial absorbing boundaries for acoustic wave equation modelling" reference.

## Code Structure

- `include/materials/helper/SochackiSpongeMaterial.h`  
  Declares the material that computes the attenuation coefficient `sochacki_damping`.

- `src/materials/helper/SochackiSpongeMaterial.C`  
  Calculates a block-wise property by:
  1. Computing normalized distances in each coordinate between the interior (`inner_*`) and outer (`outer_*`) bounds.
  2. Clamping the result to `[0,1]` and combining corner distances via `sqrt(nx*ny)` to keep the ramp smooth in corners.
  3. Evaluating the selected profile to obtain `A(x,y)` while capping at `S_max`.

- `include/kernels/core/SochackiSpongeDamping.h` and `src/kernels/core/SochackiSpongeDamping.C`  
  Provide a kernel that adds the residual contribution `2*A*rho*du/dt` (with optional `damping_scale` override) to the momentum balance. The Jacobian is zero for the explicit central-difference scheme.

## Input-Deck Integration

`development/aftershock2D/test_Sponge/dynamic_solve.i` now:

1. Defines reusable sponge parameters (`sponge_s_max`, profile controls) near the top of the file.
2. Instantiates `SochackiSpongeMaterial` on block `10`, which represents the sponge zone. The inner bounds track the elastic block extents while the outer bounds match the generated mesh limits.
3. Adds `SochackiSpongeDamping` kernels for `disp_x` and `disp_y` on block `10`, using the density material property so the damping scales correctly with local mass.
4. Removes the previous Lysmer dashpot boundary conditions; the sponge layer now handles absorption.

Run `./farms-opt -i development/aftershock2D/test_Sponge/dynamic_solve.i` to confirm the setup and tune the parameters if additional damping is needed.
