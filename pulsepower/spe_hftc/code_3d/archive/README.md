# Archived inputs (moved here 2026-10-25)

`code_3d/` keeps the pair the 3D production runs use: `elasticity_absorbing.i` (main app: Newmark elasticity, pulse
load, Lysmer dashpots on the four side faces, holder u_z = 0 on the end faces, energy and extent postprocessors, dt
schedule) and `fracture.i` (AT1 phase-field sub-app, homogeneous Gc). The files here are the 3D twins of the archived
2D variants (`../../code_2d/archive/`), none of which has been run on the production mesh.

| file | what it is | status |
| --- | --- | --- |
| `elasticity.i` | free-side-face predecessor of `elasticity_absorbing.i` (same physics without the dashpots) | smoke-mesh tests only |
| `fracture_gaussian.i` | sub-app with an element-wise Gaussian Gc field (cv 0.12, floor 0.3 Gc, seed 1) | smoke-mesh tests only |
| `fracture_gaussian_viscous.i` | the Gaussian sub-app with a rate-dependent (viscous) damage evolution | not run |
| `fracture_weibull.i` | sub-app with a Weibull Gc field (weakest link) | not run |

Rerunning one of them: MOOSE resolves relative paths (`mesh_file`, the sub-app `input_files`) from the folder of
the input file, so these do not run from `archive/`; copy the file back to `code_3d/` first.
