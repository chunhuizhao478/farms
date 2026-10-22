# Vista sweep: 2D section, Gaussian toughness scatter, 20-pulse trains

Four cases share one Grace node (`gg`, 144 cores, 36 ranks each, launched with
`ibrun -n 36 -o <offset>`), all with the strength-derived toughness (2721 J/m^2)
scattered 12 % and the quench after 25 us of each 60 us window:

| case | eta | zeta_pulse | purpose |
| --- | --- | --- | --- |
| `eta0p4_att` | 0.4 | 1e-8 s | lower end of the crack-producing range |
| `eta0p5_att` | 0.5 | 1e-8 s | the candidate baseline |
| `eta0p6_att` | 0.6 | 1e-8 s | upper end |
| `eta0p5_noatt` | 0.5 | 0 | reference without attenuation (face spall expected) |

Outputs per case in `cases_vista/<name>/`: `out.csv` (every step; `fracture_extent_core`
is the radial-crack measure, `fracture_extent` the farthest damaged node), `out.e`
(damage, particle velocity, toughness field, every 40 steps), checkpoints `out_cp/`.

Prerequisites on Vista: `$WORK/farms_cdms/farms-opt` built against `$WORK/moose-94b6528`
(see `$WORK/build_logs/`), the environment file `$WORK/build_logs/env_moose.sh`, and the
section mesh `../mesh/spe_section_h0p3.msh` (copied; `.msh` files are git-ignored).
No email notifications in this script.

## `submit_vista_att.sbatch` (the two cases requested on 2026-10-21)

| case | input | faces | eta | zeta_pulse | pulses |
| --- | --- | --- | --- | --- | --- |
| `att_eta0p5_free` | `elasticity.i` | free | 0.5 | 1e-8 s | 20 |
| `att_eta0p5_absorb` | `elasticity_absorbing.i` | Lysmer dashpots on the four side faces | 0.5 | 1e-8 s | 20 |

72 ranks each on one `gg` node. The same two cases run on Discovery from
`../cases/discovery_att_eta0p5_{free,absorb}/` (20 cores each, slower).

## Build and the first production run (2026-10-22)

Build chain on Vista: PETSc 3.25.4 on the login node; libMesh + WASP as job 1089062 on `gg`
(`--disable-netgen`, MOOSE_JOBS=64); farms as job 1089064; input check
`farms-opt -i elasticity_absorbing.i ... --check-input` = Syntax OK. The gcc14 OpenMPI libtool
archives mix two UCX versions; the workaround (a shadow lib dir without `.la` files, patched into
the PETSc link variables and reported by the wrappers through `OMPI_LDFLAGS`) is documented in
`$WORK/build_logs/README_vista_build.md` and applied by `$WORK/build_logs/env_moose.sh`.

Job 1089067 (one `gg` node, 72 ranks per case, 12 h) started 02:45 PDT. Pace: the absorbing case
covers about 120 us of the train per hour (pulse 7 after 3 h, so 20 pulses in roughly 10 h); the
free case ran 25 us per 10 min until pulse 4 and then slowed to about 20 us per hour (corner
fragments at 60-97 m/s keep the adaptive step far below `dt_max`). A job that ends by wall time is
continued with `sbatch --export=ALL,RECOVER=1 submit_vista_att.sbatch` (the session watcher does
this, up to three times). Results are tabulated in `../../local_runs/README.md`.

Status 2026-10-22 10:00 PDT: `att_eta0p5_free` was stopped by the user at 389 us (pulse 7, "not giving
correct results": corner diagonals and a fragmented sample) by killing its `prterun` on the node
(`$WORK/build_logs/kill_free_case.sh`; on Vista `ibrun` runs OpenMPI's prterun, so the two cases are not
separate SLURM steps). `att_eta0p5_absorb` continues alone in job 1089067. The sbatch now takes
`ONLY=<case>` to launch or continue a single case.

## Homogeneous-toughness twins (job 1090509, submitted 2026-10-22 13:08 PDT)

`submit_vista_homog.sbatch` runs two 20-pulse trains on one Grace node (72 ranks each), both with
the constant toughness of `fracture.i` (Gc = 2721 J/m^2, no scatter), absorbing side faces,
zeta_pulse = 1e-8 s and the quench: `att_homog_eta0p4_absorb` (eta = 0.4) and
`att_homog_eta0p5_absorb` (eta = 0.5). They are the user's test of the near-borehole picture
(2026-10-22): the local single-pulse eta = 0.4 run gave only a uniform damage ring (d = 0.32,
frozen after 6 us), so these cases ask whether repeated pulses grow that ring or nucleate
cracks, and what the homogeneous eta = 0.5 medium does against the Gaussian production case.
Same RECOVER=1 / ONLY= switches as `submit_vista_att.sbatch`; a session watcher
(`vista_watch_homog.log`) resubmits unfinished cases if the 12 h limit ends the job first.

## Volumetric-deviatoric split twin (job 1090603, submitted 2026-10-22 13:58 PDT)

`submit_vista_voldev.sbatch` runs `att_eta0p5_voldev_absorb`: the production design (Gaussian Gc,
eta = 0.5, absorbing faces, zeta_pulse = 1e-8 s, quench, 20 pulses) with
`Materials/elasticity/decomposition=VOLDEV` instead of SPECTRAL, 72 ranks on half a Grace node
(node i614-141, ~3 s per step at the start). Purpose (2026-10-22): the spectral split returns zero
driving force during the confined-compression peak at the wall and only damages through the late
hoop tension, which leaves a thin ring; the Amor split counts the shear energy of the peak (threshold
158 MPa radial stress, exceeded by the 222 MPa pulse out to 2-4 mm) and doubles the driving force in
the tensile stage, so it should give the several-mm crushed zone of the photos. Its implied unconfined
compressive strength is 73 MPa against the measured 124 MPa (ASTM D695), so the crushed band is an
upper bound. Watch item: whether the softer collar shortens or suppresses the radial cracks. A local
one-pulse twin runs on the Mac (`local_runs/section_h0p3_gauss_eta0p5_voldev_absorb_1pulse.*`).
Watcher `vista_watch_voldev.log` resubmits with RECOVER=1 if the 12 h limit ends the job first.

### Status read at 14:41-15:00 PDT, 2026-10-22 (ssh master down since 15:06 PDT, watchers paused)

- 1089067 (production): t = 1138.5 us in pulse 19 at 14:41 with 11:56 h elapsed, so the job hit
  the 12 h wall at about 14:45; 61.5 us (pulse 20) remain. Total extent 21.61 mm (capped at the side
  face since pulse 16), core extent 10.85 mm, crack energy 356 J/m, pressure work 3.97 kJ/m, dashpot
  work 10.05 kJ/m. The watcher submits the recover tail (RECOVER=1, ONLY=att_eta0p5_absorb) as soon
  as the ssh master is back.
- 1090509 (homogeneous twins) at 15:00: eta 0.4 at 296.7 us (pulse 5): extent 5.03 mm, crack energy
  85 J/m; eta 0.5 at 248.7 us (pulse 5): extent 8.37 mm, crack energy 132 J/m. The paper has 7.97 mm
  after 5 pulses and the Gaussian production case 11.04 mm. The homogeneous eta 0.4 ring (d = 0.32
  after pulse 1) does ratchet into cracks over the pulses.
- 1090603 (VOLDEV) at 14:56: t = 84.9 us (pulse 2): core extent 10.78 mm, total extent 32.3 mm
  (a corner), crack energy 611 J/m, v_max 134 m/s. Reading in `local_runs/README.md` (first look
  at 2 and 5 us): recommended to cancel.
- Watcher fix (15:40 PDT): `squeue -j` returns rc 1 ("Invalid job id specified") once a finished
  job is purged from the queue; the watchers read that as an ssh failure and retried forever, so
  the production tail was never submitted. They now treat rc 255 as the ssh failure and "invalid
  job id" as "left the queue" (then sacct decides).
