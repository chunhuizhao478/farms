# 2D section, eta = 0.5, attenuation zeta_pulse = 1e-8 s, side faces: absorb

Input: `../../elasticity_absorbing.i` + `../../fracture_gaussian.i` (Gaussian 12 % toughness scatter, seed 1),
0.3 mm mesh, l = 0.6 mm, Gc = 2721 J/m^2, 20 pulses of 60 us (25 us undamped except the
1e-8 s attenuation, then the quench). Lysmer dashpots on xmin/xmax/ymin/ymax; their work is `damping_work` in the CSV and enters `full_input_energy`.
Outputs: `out.csv` every step (core and total extents, damaged areas, energies, v_max), `out.e` every 40 steps, `out_cp/` checkpoints.
Twin case on Vista: `../../cases_vista/submit_vista_att.sbatch` (72 ranks each, same inputs).
