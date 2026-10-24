# 3D cube cases on TACC Vista

The 3D twin of the 2D section production design: homogeneous Gc, eta = 0.45, absorbing side faces.
Inputs: `../elasticity_absorbing.i` (main app; the 3D free-face `../elasticity.i` plus Lysmer dashpots
on xmin/xmax/ymin/ymax with the dashpot work in the energy budget, holder u_z = 0 on zmin/zmax kept)
and `../fracture.i` (AT1 sub-app, constant Gc = 8 l ft^2/(3E) = 2721 J/m^2 at l = 0.6 mm). Mesh
`../../mesh/spe_cube_h0p3.msh`: 14.6 M tets, 2.4 M nodes, 0.3 mm inside r <= 28 mm over the full
height, 1 mm in the corners (built locally with gmsh/HXT, 747 MB ASCII, copied to Vista on 2026-10-23;
not in git). Same loading as the 2D case (P_peak 214 MPa at the gap, 0.3 us rise, 6 us to 10 %,
zeta_pulse = 1e-8 s during the first 25 us of each 60 us window, 1e-5 s quench after it, 20 pulses),
with the 3D specifics of the plan: discharge gap centred at z = Lz/3 so the pressure decays along the
hole (139 MPa at the near end, 114 MPa at the far end for eta 0.3; the same ratios at 0.45), and the
corner node pinned in-plane.

Differences from the 2D section that matter when comparing: the 2D section applies the gap pressure
over the whole (unit) thickness, so the 3D extent at the gap plane is expected below the 2D 20.0 mm at
pulse 20; `fracture_extent` is the max over the whole cube of the radial distance of d >= 0.5 (core
measure inside r <= 15 mm as in 2D); `damping_work` now includes the z component of the dashpot force.

## Workflow

1. `submit_vista_3d_timing.sbatch` (8 nodes, 1152 ranks, 1.5 h): splits the mesh for 1152 ranks
   (`../../mesh/spe_cube_h0p3_split1152.cpr`, once) and runs the first 3 us of pulse 1 (about 50
   steps) with Exodus frames at INITIAL and FINAL. Read from `timing_8n/run.log`: seconds per step
   (steady part at dt = 1e-7), memory per node, time of one Exodus write of the full mesh, Picard
   iterations per step; from `timing_8n/out.csv`: p_wall_max 214 MPa, external_work and
   damping_work both nonzero, borehole_area 1.2813e-3 m^2.
2. Size the chain: 20 pulses are about 15 000 steps (618 per pulse without cutbacks). At the measured
   step cost pick the node count (the split must be redone for another rank count: the scripts name
   the split by the rank count and build it when missing) and the number of 24 h jobs.
3. `submit_vista_3d_eta0p45.sbatch` (8 nodes, 24 h; `RECOVER=1` continues from the newest checkpoint,
   100-step interval, 2 kept, 1152 files each). Case folder `cube_homog_eta0p45_absorb/`: `out.csv`
   every step, `out.e` every 100 steps (d and vel_mag only: 19 MB per field and frame), `run.log`
   appended across the recover chain.

Outputs are read with the 2D tools: `tools/pulse_table.py` on `out.csv` (recover-aware) gives the
extent at the end of each pulse; the Exodus frames need a 3D reader (slices at the gap height and at
the ends) instead of `exo_last_step.py`.

Local smoke test (Mac, 8 ranks, coarse mesh `spe_cube_smoke.msh` 0.15 M tets, l = 3 mm):
`PATH=$PWD/../../tools/jit_shim:$PATH mpiexec -n 8 farms-opt -i ../elasticity_absorbing.i mesh_file=../../mesh/spe_cube_smoke.msh l=3e-3 Np=1 eta=0.45 zeta_pulse=1e-8 Executioner/end_time=6e-7`
(the conda env's bin must be on PATH for the FParser JIT; result in `local_runs/README.md`).

## Timing run on 8 nodes (job 1093056, 2026-10-23, read at 13:40 PDT)

Production mesh, 14.59 M tets, 2.38 M nodes; 7.15 M displacement dofs, 2.38 M damage dofs; 1152 ranks (6.2 k displacement
dofs per rank). First 3 us of pulse 1, 49 steps, 48.5 min of wall time: 10 min for the mesh split (8 ranks, 1.7 GB), 45 s
of setup, 72 s per full-mesh Exodus frame (serialised to rank 0), 38 min of stepping.

| phase (step, dt) | fixed-point its | main app Newton / GMRES | main solve | step wall |
| --- | --- | --- | --- | --- |
| rise, dt 30 ns (steps 6-27) | 1 | 2 / 19 | 21 s | about 30 s |
| plateau before cracking, dt 100 ns (28-35) | 1 | 1 / 5 | 13 s | about 20 s |
| crack nucleation at the wall, dt 100 ns (36-49) | 4 | 4-5 / 15-20 (summed) | 34-44 s | about 60 s |

The sub-app needs 0 Newton iterations until the wall reaches the strength (step 34, t = 1.5 us) and 7-9 per step
(summed over the four fixed-point iterations) once cracks nucleate; its solve time is not in the live print (sub-apps do
not live-print), the step wall times above are the run total minus the printed main-app solves. Picard converges in one
iteration without damage evolution and in four with it, as in 2D. Pulse law reproduced at the wall (p_wall_max follows the
double exponential to within 0.1 %), damping work zero before the wave reaches the faces, crack area 31 mm2 at 3 us.

Projection with the 2D step count (12.5 k steps for 20 pulses, dt_max 1e-7): roughly 250 cracking steps and 350 quiet
steps per pulse give about 7 h per pulse, 140 h for 20 pulses on 8 nodes (about 1100 node-hours of the 25 125 SU balance),
i.e. three 48 h segments (qgg allows 2 days, 32 nodes per job). Sizing runs on 4 nodes (job 1093797) and 16 nodes
(job 1093798) decide the node count; both print the perf graph at step 25 and at the end (the inputs now carry a
PerfGraphOutput every 500 steps, level 2, for the main and the sub app, so that production segments killed by the
wall limit still leave a cost profile). The production script was switched to 48 h segments.

## Timing run on 4 nodes (job 1093797, 576 ranks, 2026-10-23 17:30 PDT) and the cost breakdown

Same 49 steps to 3 us: 3196 s of stepping against 2078 s on 8 nodes, i.e. 8 nodes are 1.54 times faster than 4 (77 %
parallel efficiency). Level-2 perf graph of the main app over the 63 fixed-point iterations (4 nodes):

| section | share | note |
| --- | --- | --- |
| nonlinear solve (AD assembly dominates: Jacobian 7.6 s, residual 4.0 s per call; linear solve small) | 58 % | 1.8 Newton its per fixed-point it |
| fixed-point residual-norm checks (two extra residual evaluations per iteration) | 14 % | avoidable with a custom_pp criterion |
| postprocessors (energies, extents, areas over the whole mesh every step) | 9 % | |
| fracture sub-app (0 to 9 Newton its here; grows with crack growth) | 6 % | |
| transfers, console, csv | 4 % | |
| Exodus: 133 s for the first frame (mesh), 9 s per later frame | 1 % at one frame per 100 steps | |

Memory 18.5 GB on rank 0 (Exodus serialisation), far below the node limit. The sub-app perf graph is summed over ranks
and not usable for timing. Decision pending the 16-node point (job 1093798): the planned production is 32 nodes (the
per-job maximum) with the ring-down cap dt_quiet = 5e-7, which should fit the 20 pulses in a single 48 h job
(about 35 h, 1100 node-hours) without a recover chain; dt_max = 2e-7 in the growth phase and the custom_pp fixed-point
criterion are the further levers if a day is required, both to be checked in 2D first (submit_vista_dtcheck.sbatch).

## Timing run on 16 nodes (job 1093798, 2304 ranks, finished 2026-10-24 00:05 CDT) and the node count

Same 49 steps to 3 us, 1614 s of stepping (`Transient::PicardSolve`, 63 fixed-point iterations). Scaling of the
stepping time:

| nodes | ranks | stepping | s per step | speed-up vs 4 nodes | parallel efficiency |
| --- | --- | --- | --- | --- | --- |
| 4 | 576 | 3196 s | 65 | 1 | |
| 8 | 1152 | 2078 s | 42 | 1.54 | 77 % |
| 16 | 2304 | 1614 s | 33 | 1.98 | 50 % (64 % from 8 to 16) |

Per call on 16 nodes: Jacobian 3.2 s (7.6 s on 4), residual 2.4 s (4.0 s on 4): the residual assembly scales worse
than the Jacobian, so the fixed-point norm checks grow to 17 % of the stepping. First Exodus frame 170 s (133 s on 4
nodes; the serialisation to rank 0 grows with the rank count), console 1 s per step, 18.7 GB on rank 0. The queue wait
of this 16-node, 1.5 h job was 8 h (submitted 15:23, started 23:2x CDT). Extrapolating the efficiency decay, 32 nodes
would gain 10-20 % over 16 at twice the node-hours and a longer wait.

Decision: 16 nodes (`spe_cube_h0p3_split2304.cpr` exists, 2.0 GB; the production script defaults to -N 16 -n 2304).
Projection with the 2D step counts for 20 pulses (240 rise steps, 5129 in the growth phase 0.3-25 us of which 3707
advance the crack, 7136 in the ring-down 25-60 us) and the 16-node costs (33 s per quiet step, about 46 s per step
with crack growth, the 8-node ratio 60/42):

| schedule | steps | 16 nodes | 32 nodes (assumed 1.15x) | segments of 48 h |
| --- | --- | --- | --- | --- |
| as is (dt_quiet = dt_max = 1e-7) | 12 500 | about 130 h | 110 h | 3 |
| dt_quiet = 5e-7 (ring-down 1400 steps) | 6 800 | about 75 h | 65 h | 2 (one recover) |
| dt_quiet = 5e-7 and dt_max = 2e-7 | 4 200 | about 45 h | 39 h | 1 (marginal on 16 nodes) |

This supersedes the 35 h estimate above, which extended the 77 % efficiency to 32 nodes. Under 24 h is out of reach
at this mesh size: 32 nodes with both caps and the custom_pp fixed-point criterion (removes the 17 % norm checks)
would land near 33 h. dt_quiet = 5e-7 costs nothing physically (kinetic energy below 3 % of the window peak at 25 us);
dt_max = 2e-7 in the growth phase changes the crack increments per step and must pass the 2D check first
(`code_2d/cases_vista/submit_vista_dtcheck.sbatch`: extents at pulses 5/11/15/20 and crack energy against the
baseline, plus the step count and wall time, which also tell whether larger steps cost more Picard iterations).
Recommended path: 16 nodes with dt_quiet = 5e-7, two segments through the watcher; switch to one segment only if the
2D check passes for dt_max = 2e-7.

Submission (from this folder, on Vista): `sbatch --export=ALL,EXTRA="dt_quiet=5e-7" submit_vista_3d_eta0p45.sbatch`;
continuation `sbatch --export=ALL,RECOVER=1,EXTRA="dt_quiet=5e-7" submit_vista_3d_eta0p45.sbatch` (the watcher
`vista_watch_3d.sh JOBID NNODES` passes its own EXTRA through). EXTRA is appended to the command line, so any
input parameter can be overridden there; the input defaults are unchanged (dt_quiet = dt_max).
