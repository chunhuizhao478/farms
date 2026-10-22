# Discovery run: 2D section, Gaussian toughness scatter, eta = 1, 20-pulse train

What it computes: the plane-strain section of the acrylic cube at the gap height
(`../../elasticity.i` + `../../fracture_gaussian.i`), 0.3 mm mesh, phase-field
length 0.6 mm, toughness 2721 J/m^2 (from ft = 69 MPa) scattered element-wise with
a truncated Gaussian of 12 % standard deviation (bounds 0.3x and 2x the mean, seed 1),
all 144 J per pulse reaching the wall (eta = 1), twenty 60 us windows (pulse +
25 us free + 35 us quench), i.e. end_time = 1.2 ms, about 15 000 steps.

Outputs (in this folder, file base `out`):
- `out.e`: every 40th step and the last one, nodal `d` (damage), `vel_x`, `vel_y`,
  `vel_mag` (particle velocity magnitude), elemental `Gc_var` (the random toughness
  realization) and `pulse_load_aux` (applied pressure field).
- `out.csv`: every step: energies (input, kinetic, elastic, crack), `fracture_extent`
  (d >= 0.5) and `fracture_extent_090`, `damaged_volume` (area per unit thickness here),
  `crack_area` (crack energy / Gc, a length per unit thickness here), `p_wall_max`,
  `v_max` (maximum particle velocity), `borehole_area` (2*pi*rb = 2.6075e-2 m).
- `out_cp/`: checkpoints every 100 steps for `--recover`.

Submitting (the script is filled for user zhaochun: account scec_608, modules
usc/13.3.0 gcc/13.3.0 openmpi/5.0.5 openblas/0.3.28 python/3.11.9, farms-opt at
/home1/zhaochun/projects/farms_cdms/farms-opt, launched with mpirun inside the allocation):
1. Copy `pulsepower/spe_hftc/{code_2d,mesh/spe_section_h0p3.msh}` to the Discovery clone
   `/home1/zhaochun/projects/farms_cdms/pulsepower/spe_hftc/` (`.msh` files are git-ignored; 3 MB).
2. `sbatch submit_discovery.sbatch`; if it times out, `sbatch --export=ALL,RECOVER=1 submit_discovery.sbatch`.
First submission: job 13293853 on 2026-10-21, cancelled by the user after two pulses (1 h 50 min): at eta = 1 the section spalls at the free faces during the first pulse (see local_runs/discovery_snapshot_damage.png). Resubmitted as discovery_gaussian_eta0p3 and discovery_gaussian_eta0p5.

Expected cost: the Mac does 1.5 s per step on 8 ranks for this mesh; 32 Discovery
cores should be around 0.4 s per step, so roughly 2 hours for the train.

Variants from the same folder (append to ARGS or copy the folder):
- Weibull instead of Gaussian: `MultiApps/fracture/input_files=fracture_weibull.i`
- uniform toughness: drop the `MultiApps/...` term
- other energy fraction: `eta=0.5`
- no quench (energy-budget check): `zeta_quench=0`
