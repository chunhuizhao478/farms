# Damage-paper planar fault, S = 0.2, adaptive mesh (Zhao et al. 2024)

Self-contained reproduction of the continuum damage-breakage (CDBM) planar-fault case
from Zhao et al. 2024. Everything needed to rerun it is in this folder.

## The problem

A 2D plane-strain dynamic rupture on a planar fault, coupled to the CDBM bulk rheology.
Elastodynamics and the damage/breakage evolution are solved as two MOOSE apps, coupled
as a `TransientMultiApp`: `test_planarfault_main.i` drives, `test_planarfault_sub.i`
integrates the damage ODEs and copies alpha and B back every timestep.

| | |
| --- | --- |
| domain | 20 km x 20 km, x and y in [-10 km, +10 km] |
| fault | planar, along y = 0, full width |
| friction | linear slip-weakening, mu_s spatially set, mu_d = 0.1, Dc = 0.4 m |
| background stress | sigma_xx = -135 MPa, sigma_yy = -120 MPa, sigma_xy = 70 MPa |
| strength ratio | S = 0.2 |
| elastic moduli | lambda_o = shear_modulus_o = 3.204e10 Pa, density 2670 kg/m^3 |
| time integration | central difference (main, lumped), SSP-RK3 order 3 (sub) |
| timestep | dt = 5e-4 s, fixed |
| duration | t = 3.0 s, i.e. 6000 steps |
| output | Exodus every 200 steps -> 31 frames at t = 0.0, 0.1, ... 3.0 |

### Mesh

`mesh/planarfault2_uniform_adaptive.msh`: 373,068 nodes, 746,940 elements. Triangles,
25 m within 3 km of the fault, grading to 125 m at the domain edge (gmsh `Distance` +
`Threshold` fields; see the `.geo`).

MOOSE splits the mesh along the fault with `BreakMeshByBlockGenerator`, so the Exodus
reports slightly different counts — 373,869 nodes and 744,882 triangles — because the
801 fault nodes are duplicated and the mesh file's 1D line elements are dropped. That
difference is expected, not a sign of a different mesh.

**The `.msh` is committed, not just the `.geo`.** Other examples here ship only a `.geo`.
This one ships both deliberately: gmsh output is version-sensitive, and a regenerated
mesh would have different node counts, which would break the element-count checks in
`verify_fingerprint.py`. The committed mesh was produced with gmsh 4.1 (file format 4.1).

## Build

This example needs `ComputeDamageBreakageStressv2`, `FDCompVarRate2` and
`BreakageVarForcingFuncDevOld`, all present on this branch. A normal app build suffices:

```bash
git clone --recurse-submodules https://github.com/chunhuizhao478/farms.git
cd farms
export CC=mpicc CXX=mpicxx FC=mpif90 F90=mpif90 F77=mpif77
export METHOD=opt METHODS=opt MOOSE_JOBS=16
make -j 16
```

Build MOOSE first if you do not already have one; `moose` is a submodule pinned at
`8fd96527` (2024-07-02). See `SYNTAX_MIGRATION.md` for why that pin matters.

## Run

```bash
cd examples/damage_paper_zhaoetal2024
mpirun -np 128 ../../farms-opt -i test_planarfault_main.i
```

`run_expanse.sb` is a worked SLURM script with preflight assertions that fail *before*
consuming a node-day. Adapt `--account` and the module lines for your site.

**Cost, measured on SDSC Expanse** (2 x AMD EPYC 7742, 128 cores/node):

| | |
| --- | --- |
| wall time | 22 h 09 m on 1 node / 128 ranks |
| charge | ~2,836 core-hours |
| rate | ~10.5 s per timestep |
| output | `test_planarfault_main_out.e` ~12 GB, plus ~2 GB of snapshot files |

Budget 24 h. The reference run finished with under 2 h to spare.

## Verify

**Do not judge success by "it ran", and do not judge it by output file size.**

```bash
python3 verify_fingerprint.py test_planarfault_main_out.e
```

This compares per-frame alpha and B statistics against `reference/fingerprint_alpha_B.csv`,
measured from the original 2023 output. Element counts must match **exactly** — they are
integers set by threshold crossings, so any drift in the fields moves them. `max` and
`mean` must agree to 6 significant figures.

The two signatures that matter most:

- **alpha reaches exactly 1.0** from t = 1.7 s onward.
- **B is identically zero through t = 1.3 s**, then appears at t = 1.4 s in exactly
  **15 elements** with max **0.605578**, and ends at t = 3.0 s with **48,643** broken
  elements.

## Known pitfalls

These have each cost real time. Read them before debugging anything.

### 1. The alpha_cr calibration — the one that fails silently

`alpha_cr` (critical damage: where solid converts to granular) is a closed form in xi
whose coefficients bake in `lambda_o` and `shear_modulus_o`. A given polynomial is only
valid for the moduli it was derived for. This case runs `3.204e10`; the kernel's default
is derived for `2.73e9`.

Both `BreakageVarForcingFuncDevOld` blocks in `test_planarfault_sub.i` therefore set:

```
alphacr_calibration = stiff_3p204e10
```

If that is missing, **the run does not fail**. It produces plausible output in which
alpha caps at **0.952169427** instead of reaching 1.0, and 45,179 elements are broken at
t = 3.0 s instead of 48,643 — a 7% error with no warning. `verify_fingerprint.py` detects
this and names it.

The kernel also cross-checks the choice against `lambda_o`/`shear_modulus_o` from
`[GlobalParams]` and errors on a mismatch beyond 1%, so a wrong pairing is now loud
rather than silent.

### 2. Matching output file size proves nothing

A rerun with a differently-calibrated build produced an Exodus of 11,940,877,188 bytes
against the reference's 11,940,860,476 — **17 KB apart in 11.94 GB, 0.00014%** — while
containing materially different damage. File size tracks mesh, variable set and frame
count. It says nothing about the values inside. Use `verify_fingerprint.py`.

### 3. Deprecation warnings are benign

The run prints deprecation warnings for `Modules/TensorMechanics/Master`,
`Modules/TensorMechanics/CohesiveZoneMaster` and `Outputs/interval`. All three are still
supported at the pinned MOOSE. See `SYNTAX_MIGRATION.md`, which also gives the exact
replacements should you bump the submodule.

### 4. Never read the output while the job is still writing

On a parallel filesystem (Lustre), reading the in-progress Exodus from a login node
forces write-lock revocation and can stall the writer. Doing this once cost ~20 minutes
of a 22-hour run. Wait for the job, or copy the file first.

### 5. Reading a 12 GB Exodus needs care

An `mmap` of the output is refused under a typical login-node address-space limit, and
`scipy.io.netcdf_file` *without* mmap tries to load the whole record array into RAM.
`verify_fingerprint.py` reads frame by frame; run it as a small batch job if your login
node is tight.

## Changes from the published deck

The decks are `farms` commit `99b41e9c` (2023-11-20) verbatim, with exactly two edits:

1. **Mesh path** repointed to the sibling `mesh/` folder.
2. **`alphacr_calibration = stiff_3p204e10`** added to both breakage kernels. On the
   original branch this was the hardcoded behaviour, so this line preserves the published
   physics rather than changing it.

No syntax migration was needed — see `SYNTAX_MIGRATION.md`.

Note in particular what was **not** changed: `num_steps` stays commented out and
`Outputs/interval` stays at 200. A later commit on the `cdbm` branch flipped these to
`num_steps = 10` and `interval = 1` as a debug throttle; running that version stops after
5 ms of a 3 s rupture and writes every timestep. `run_expanse.sb` asserts both.

## Files

| file | what it is |
| --- | --- |
| `test_planarfault_main.i` | elastodynamics + slip-weakening fault, drives the MultiApp |
| `test_planarfault_sub.i` | damage/breakage evolution (SSP-RK3) |
| `mesh/*.geo`, `mesh/*.msh` | gmsh source and the committed mesh |
| `run_expanse.sb` | worked SLURM script with preflight assertions |
| `verify_fingerprint.py` | pass/fail check against the published fields |
| `reference/` | the fingerprint and its provenance |
| `SYNTAX_MIGRATION.md` | MOOSE deprecation status at the pinned submodule |
| `tests/` | standalone equivalence test for the alpha_cr calibration |
