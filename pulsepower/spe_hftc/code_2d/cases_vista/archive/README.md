# Archived job scripts (moved here 2026-10-25)

`cases_vista/` keeps one job script, `submit_vista_eta0p45.sbatch`: the homogeneous-Gc eta 0.45 case
(`att_homog_eta0p45_absorb`), the calibrated 2D reference and the twin of the 3D production case. The scripts
below are the records of the other Vista runs; the case folders and the results stay in `cases_vista/` and in its
README. All of them `cd` to `code_2d` and write to `cases_vista/<case>/`, so they can be submitted from here
(`sbatch archive/<script>` from `cases_vista`), except that the inputs `elasticity.i` and `fracture_gaussian.i` they
name now live in `code_2d/archive/` and must be copied back to `code_2d/` first.

| script | cases | inputs | jobs |
| --- | --- | --- | --- |
| `submit_vista_att.sbatch` | `att_eta0p5_free` (free side faces), `att_eta0p5_absorb` (Gaussian Gc, the first production candidate) | `elasticity.i` / `elasticity_absorbing.i` + `fracture_gaussian.i` | 1089067, continued as 1092063 |
| `submit_vista_homog.sbatch` | `att_homog_eta0p4_absorb`, `att_homog_eta0p5_absorb` | `elasticity_absorbing.i` + `fracture.i` | 1090509, continued as 1092069 |
| `submit_vista_voldev.sbatch` | `att_eta0p5_voldev_absorb` (volumetric-deviatoric split test) | `elasticity_absorbing.i` + `fracture_gaussian.i` | 1090603 |
| `submit_vista_eta0p425.sbatch` | `att_homog_eta0p425_absorb` (threshold check) | `elasticity_absorbing.i` + `fracture.i` | 1093037 |
| `submit_vista_sweep.sbatch` | early free-face sweep, 36 ranks per case | `elasticity.i` + `fracture_gaussian.i` | superseded by the absorbing cases |
| `submit_vista_dtcheck.sbatch` | `att_homog_eta0p45_absorb_dtq5`, `..._dtq5_dtm2` (step-size check for the 3D schedule) | `elasticity_absorbing.i` + `fracture.i` | not submitted; runs as is from here |

The Gaussian eta 0.45 case (`att_eta0p45_absorb`, job 1092158) was launched by the kept script, whose Gaussian run
line is now commented out.
