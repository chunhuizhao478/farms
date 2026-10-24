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

### Vista system event on 2026-10-22 (17:00-18:40 CDT) and the resubmissions (21:05 PDT)

- All three runs of jobs 1090509 (homogeneous twins, node i618-022) and 1090603 (VOLDEV, node
  i614-141) died at the same second, 18:39:32 CDT, each with one rank killed by signal 7 (Bus error)
  while writing output or computing a residual. Their CSV, Exodus and checkpoint files stop at
  17:12-17:29 CDT, i.e. the solvers hung for over an hour before the kill, and the ssh master to the
  login node had dropped at 17:06 CDT. Two nodes failing in the same second on the shared $WORK
  file system reads as a Lustre event, not a model problem. Slurm reports both jobs as COMPLETED
  (the job scripts exited normally after ibrun died).
- Final states: homogeneous eta 0.4 at 317.3 us (pulse 6): extent 6.69 mm, crack energy 105 J/m;
  homogeneous eta 0.5 at 270.5 us (pulse 5): extent 10.16 mm, crack energy 150 J/m; VOLDEV at
  89.5 us (pulse 2): core extent 10.78 mm, total extent 32.3 mm, crack energy 615 J/m.
- 1089067 (production) was unaffected: it timed out at 16:45 CDT at 1143.9 us (pulse 20, 3.9 us
  into the last pulse). The watcher submitted the recover tail as job 1092063 at 22:59 CDT
  (checkpoint 12100 of 16:45 CDT, intact); it is stepping from 1145 us to the 1200 us end on
  i615-014, about 35 minutes of work.
- Homogeneous twins resubmitted with RECOVER=1 as job 1092065 (23:01 CDT): eta 0.4 recovered from
  checkpoint 3200 (17:12 CDT, before the hang), but eta 0.5 aborted at once because checkpoint
  2800 (17:29 CDT, the onset of the hang) holds 5 of its 73 files. Fix: 1092065 cancelled, the
  partial `out_cp/2800-*` removed so that --recover picks the complete checkpoint 2700, both twins
  resubmitted together as job 1092069 (23:05 CDT, node i617-141). Rule for next time: after a crash,
  count the files of the newest checkpoint (73 per checkpoint at 72 ranks) before recovering.
- VOLDEV not resubmitted: the one-pulse result (local_runs/README.md) and the 89.5 us state rule it
  out as a production candidate; its checkpoints (0800-1000) stay on Vista if a continuation is
  wanted.
- Watchers: `vista_watch.sh 1092063` and `vista_watch_homog.sh 1092069`; the VOLDEV watcher is off.

## eta = 0.45 pair (job 1092158, submitted 2026-10-22 21:31 PDT, node i614-061)

`submit_vista_eta0p45.sbatch`: `att_eta0p45_absorb` (Gaussian Gc, the production design at
eta 0.45) and `att_homog_eta0p45_absorb` (homogeneous Gc), 72 ranks each on one node, 12 h, no mail.
Pulse pressure 214 MPa (222.5 at eta 0.5, 206 at eta 0.4). Watcher `vista_watch_eta0p45.sh`,
log `vista_watch_eta0p45.log`.

Why: in the extent-per-pulse chart the production case grows at 0.88 mm per pulse from pulse 5 on
(experiment 0.83 averaged over pulses 5-20), so the 3 mm overshoot is built in pulses 1-4. Lowering
eta mainly delays nucleation (first cracks in pulse 1 / 2 / 3 for Gaussian 0.5, homogeneous 0.5,
homogeneous 0.4) and trims the increments by about 10 % per 0.05 of eta (threshold-like, pressure
above about 150 MPa). Projection for eta 0.45 from the production increments scaled by 0.9:
homogeneous 7-8 / 13.5-14.5 / 16-17 / 19.5-20.5 mm at pulses 5 / 11 / 15 / 20 (measured 7.97 /
13.90 / 15.66 / 20.45); Gaussian 9.5-10 / 16-16.5 / 18.5-19 / at the face. If the Gaussian case
overshoots as projected, a Gaussian eta 0.4 is the next candidate. The experiment's flat stretch
between pulses 11 and 15 (0.44 mm per pulse) is not reachable with eta alone.

VOLDEV dropped (user decision, 21:25 PDT): neither the Vista case nor the local run is continued.

## Twins complete; eta 0.45 pair read at 08:15 PDT, 2026-10-23

Job 1092069 (homogeneous twins, recovered from checkpoints 3200 / 2700 of job 1090509) finished at
07:58 PDT after 10 h 52 min, both cases at 1200 us. Job 1092158 (eta 0.45 pair): the homogeneous
case reached 1200 us at 07:40 PDT (10.2 h), the Gaussian case at 09:27 PDT, five minutes inside
the 12 h wall (11 h 55 min), so neither needed a recover; final frame and CSV fetched at 09:36 PDT. Figures (local_runs/):
`vista_snapshots_4cases.png` is now a 2 x 3 grid (rows Gaussian / homogeneous Gc, columns eta
0.5 / 0.45 / 0.4; the empty slot holds the table below) and `vista_extent_vs_pulse_allcases.png`
carries the six cases (colour = eta, solid lines and circles = Gaussian, dashed and triangles =
homogeneous; VOLDEV removed). Per-pulse numbers from `tools/pulse_table.py`, which drops the
orphaned rows after a recover checkpoint and sums the crack-energy column across the resets.

Crack extent from the wall at the end of pulses 5 / 11 / 15 / 20 in mm (measured 7.97 / 13.90 /
15.66 / 20.45); energies per unit thickness at 1200 us:

| case | 5 | 11 | 15 | 20 | cracks > 2 mm, longest | crack energy | pressure work | dashpot work |
|---|---|---|---|---|---|---|---|---|
| Gaussian eta 0.5 (production) | 11.04 | 17.16 | 19.84 | 21.61 (side face) | 13, 21.6 mm | 388 J/m | 14.9 kJ/m | 10.6 kJ/m |
| Gaussian eta 0.45 | 13.19 | 18.09 | 19.67 | 22.31 (side face) | 10, 22.3 mm | 313 J/m | 13.6 kJ/m | 9.7 kJ/m |
| homogeneous eta 0.5 | 10.16 | 18.08 | 19.85 | 23.31 | 9, 23.3 mm | 490 J/m | 15.8 kJ/m | 10.7 kJ/m |
| homogeneous eta 0.45 | 7.56 | 16.11 | 18.95 | 20.01 | 13, 20.0 mm | 281 J/m | 13.5 kJ/m | 10.0 kJ/m |
| homogeneous eta 0.4 | 5.03 | 10.09 | 10.39 | 10.58 | 8, 10.6 mm | 158 J/m | 12.3 kJ/m | 9.6 kJ/m |

Increments per pulse in mm, pulses 2-20: homogeneous 0.45: 0.5, 2.6, 2.3, 2.2, 1.8, 1.7, 1.4, 1.3,
1.2, 1.2, 1.0, 0.8, 0.6, 0.4, 0.3, 0.2, 0.3, 0.1, 0.2; homogeneous 0.4: 0.0, 0.5, 2.4, 2.1, 1.7,
1.4, 1.0, 0.5, 0.3, 0.3, 0.2, then 0.0-0.2 (arrest); homogeneous 0.5: 1.6, 4.1, 2.7, 1.8, 1.3, 1.5,
1.5, 1.5, 1.4, 0.7, 0.7, 0.4, 0.3, 0.3, 0.6, 0.6, 0.8, 0.6, 0.9; Gaussian 0.45: 4.0, 2.8, 2.5,
1.8, 1.2, 0.8, 0.8, 0.8, 0.8, 0.6, 0.5, 0.4, 0.4, 0.3, 0.3, 0.3, 0.3, 0.3, 1.5 (the crack at 344 deg
reaches the right face in pulse 20, as the production's did in pulse 16).

Reading:

- Homogeneous eta 0.45 is the closest history: within 0.5 mm of the measurement at pulses 5 and
  20, 2.2 and 3.3 mm high at pulses 11 and 15. The model grows at a steady 1.2-1.4 mm per pulse
  from pulse 5 to 12 and slows only as the leading cracks approach the faces, whereas the
  experiment has its slow stretch between pulses 11 and 15 (0.44 mm per pulse) and then 1.0 mm
  per pulse to pulse 20.
- Homogeneous eta 0.4 arrests: the six-armed star stops at 10.6 mm and nothing moves after pulse
  12. Between eta 0.4 and 0.45 the extent at pulse 20 doubles, so eta acts as a threshold rather
  than a slope; an eta 0.425 would arrest somewhere between 11 and 16 mm and could not reach 20 mm
  at pulse 20.
- Homogeneous eta 0.5 overshoots at every reading and is still growing at 0.6-0.9 mm per pulse at
  the end (longest crack 23.3 mm at 60 deg, 2 mm short of the top face).
- The Gaussian eta 0.45 case runs ahead of the Gaussian eta 0.5 case in pulses 2-5 (13.2 vs
  11.0 mm at pulse 5) with 4 % less pressure and less crack energy (160 vs 168 J/m at pulse 5, 285
  vs 330 at pulse 18): same seed, but fewer cracks take off (10 vs 13 longer than 2 mm) and the
  crack at 346 deg gets a larger share. The maximum extent of one Gaussian realization is
  therefore a noisy measure, with a realization scatter of about 2 mm, as large as the eta effect;
  the crack energy is the robust one (homogeneous 490 / 281 / 158 J/m for eta 0.5 / 0.45 / 0.4,
  about -40 % per 0.05 of eta).
- Crack energy is 1.3-3.1 % of the pressure work in every case; the dashpots take 9.2-10.7 kJ/m.

Recommendation: calibrate eta on the homogeneous model, where eta 0.45 is the answer. If the
Gaussian design stays in the paper, run it with two more seeds at eta 0.45 and report the band
instead of moving eta against one realization. The Gaussian eta 0.4 candidate is withdrawn.

## eta = 0.425 homogeneous (job 1093037, submitted 2026-10-23 09:40 PDT)

`submit_vista_eta0p425.sbatch`: `att_homog_eta0p425_absorb` alone, 72 ranks (the same decomposition as
the eta 0.4 / 0.45 / 0.5 cases, so that the symmetry breaking of the crack star is seeded the same
way), 12 h, no mail; pulse pressure about 210 MPa. Watcher `vista_watch_eta0p425.sh`, log
`vista_watch_eta0p425.log`; it resubmits with `--recover` on a wall-time end. The second half of the
node is idle (one case per node was the user's choice).

Why: the three homogeneous runs fall into two regimes, not on a slope. eta 0.4 keeps a symmetric
six-armed star (arms 7.0-10.6 mm) that arrests after pulse 12; eta 0.45 localizes, two diagonal arms
run to 20.0 and 19.4 mm while the other four stay at 7.8-11.2 mm; the peak pressures differ by 4 %.
eta 0.425 sits inside that gap. Near such a threshold the star should stall for several pulses before
one arm breaks away, i.e. a fast-slow-fast history, which is the shape of the measurement (0.44 mm per
pulse between pulses 11 and 15, 1.0 mm per pulse afterwards). Outcomes: a stall at 11-13 mm followed
by a late localization gives a calibrated eta with a mechanism behind the measured shape; an arrest
puts the threshold between 0.425 and 0.45 and leaves eta 0.45 as the calibration, with the overshoot
at pulses 11 and 15 as the price of being above threshold. Plain interpolation (about 6 / 13 / 14.5 /
15 mm at pulses 5 / 11 / 15 / 20) is the least likely outcome. Caveats: near the threshold the mesh
asymmetry decides which arm wins and when, so one run gives the regime, not a precise curve; the 2D
section cannot carry the 3D crack front.

Result (finished 2026-10-23 20:44 PDT, 8:28 h on 72 ranks, 12 381 steps, no recover): arrest, the
eta 0.4 regime. The star has the same arms at the same angles as eta 0.4 (two leading arms at 54 and
310 deg, four at 117-246 deg, so the mesh seeds the pattern); the longest arm is 10.5 mm at pulse 10
and 11.2 mm at pulse 20 with nothing gained after pulse 16, against 10.6 mm for eta 0.4 and 20.0 mm
for eta 0.45. Crack energy 169 J/m (158 / 281 J/m for eta 0.4 / 0.45): the extra 2 % of pressure
feeds the same arms rather than a new one. The early pulses are faster than eta 0.4 (6.7 vs 5.0 mm at
pulse 5, 1.2 mm under the measurement) and the arrest comes one pulse later; there is no stall
followed by a late localization. The localization threshold therefore lies between eta 0.425 and
0.45, eta 0.45 stays the calibration, and the measured slow stretch between pulses 11 and 15 is not
something a different eta reproduces. Figures: `local_runs/vista_extent_vs_pulse_homog_with_eta0p425.png`
(`tools/plot_extent_vs_pulse_homog.py STEM --with-0p425`, working chart; the paper figure keeps the
three calibration values) and `local_runs/vista_snapshots_eta0p4_0p425_0p45_pulse20.png`
(`tools/fig_snapshot_row.py`). Window-end frames (pulses 4, 8, 10, 12, 16, 20) were extracted on
Vista with `$WORK/tmp_npz/exo_window_frames.py`.
