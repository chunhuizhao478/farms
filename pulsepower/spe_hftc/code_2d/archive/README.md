# Archived inputs (moved here 2026-10-25)

`code_2d/` keeps the pair every current case runs: `elasticity_absorbing.i` (main app: Newmark elasticity, pulse
load, Lysmer dashpots on the side faces, energy and extent postprocessors, dt schedule) and `fracture.i` (AT1
phase-field sub-app, homogeneous Gc). The files here are the earlier variants, kept for the record of the runs that
used them.

| file | what it is | runs that used it |
| --- | --- | --- |
| `elasticity.i` | free-side-face predecessor of `elasticity_absorbing.i` (same physics without the dashpots) | Discovery cases in `../cases/discovery_*`, Vista `att_eta0p5_free`, the Vista sweep |
| `fracture_gaussian.i` | sub-app with an element-wise Gaussian Gc field (cv 0.12, floor 0.3 Gc, seed 1) | Vista `att_eta0p5_absorb` (first production candidate), `att_eta0p45_absorb`, the VOLDEV test, the sweep; Discovery Gaussian cases |
| `fracture_gaussian_viscous.i` | the Gaussian sub-app with a rate-dependent (viscous) damage evolution | local tests only |
| `fracture_weibull.i` | sub-app with a Weibull Gc field (weakest link) | local tests only |

Rerunning one of them: MOOSE resolves relative paths (`mesh_file`, the sub-app `input_files`) from the folder of
the input file, so these do not run from `archive/`. Copy the file back to `code_2d/` and launch as before, e.g.
`MultiApps/fracture/input_files=fracture_gaussian.i` on the command line of `elasticity_absorbing.i`. The job
scripts in `../cases_vista` and `../cases` that used these files still name the old location.
