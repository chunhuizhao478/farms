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
