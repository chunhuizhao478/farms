# 2D section, eta = 0.5, attenuation zeta_pulse = 1e-8 s, side faces: free

Input: `../../elasticity.i` + `../../fracture_gaussian.i` (Gaussian 12 % toughness scatter, seed 1),
0.3 mm mesh, l = 0.6 mm, Gc = 2721 J/m^2, 20 pulses of 60 us (25 us undamped except the
1e-8 s attenuation, then the quench). Free side faces (the physical specimen); the comparison case has absorbing faces.
Outputs: `out.csv` every step (core and total extents, damaged areas, energies, v_max), `out.e` every 40 steps, `out_cp/` checkpoints.
Twin case on Vista: `../../cases_vista/submit_vista_att.sbatch` (72 ranks each, same inputs).

Snapshot after two pulses (t = 108 us, 2026-10-22): pulse 1 gives 7 radial cracks to 2.65 mm from the wall with sub-critical corner streaks (d <= 0.24); pulse 2 extends the hole cracks to 6.48 mm by 76 us and then, 20 us after its onset, breaks the four corner diagonals (d = 1 along |x| = |y| at 15-16 mm, about 10 mm from the faces; no face spall). Crack energy 207 J/m, of which ~90 J/m are the corner cracks. Picture: `local_runs/discovery_att_eta0p5_free_vs_absorb_pulse2.png`.

Status 2026-10-22 09:53 PDT: cancelled by the user after 10 h 55 m (too slow at 20 cores; the Vista
twin `cases_vista/att_eta0p5_free` runs the same input three times faster). Partial outputs stay in
this folder on Discovery.
