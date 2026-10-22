# Local 2D section runs (Mac, 8 ranks, JIT shim on PATH)

Section mesh `spe_section_h0p3` (l = 0.6 mm, Gc = 2721 J/m^2 from ft = 69 MPa),
one pulse, first 10 us only (no quench yet). Command:

    PATH=../tools/jit_shim:$PATH mpiexec -n 8 farms-opt -i elasticity.i Np=1 Executioner/end_time=1e-5 Outputs/file_base=<name> [eta=1.0 | Gc_const=350 --allow-unused | MultiApps/fracture/input_files=fracture_weibull.i | MultiApps/fracture/input_files=fracture_gaussian.i]

| run | eta | Gc (J/m^2) | result |
| --- | --- | --- | --- |
| `section_h0p3_eta0p3_10us` | 0.3 | 2721 | thin uniform damage ring, d_max = 0.18 within 0.55 mm of the wall; crack energy 5.6 J/m; no node above d = 0.5 |
| `section_h0p3_eta1_10us` | 1.0 | 2721 | crushed collar: d_max = 0.97, d > 0.5 out to 1.76 mm from the wall, d > 0.9 to 0.75 mm; crack energy 77 J/m; two incipient radial cracks (1-1.5 mm) break the symmetry |
| `section_h0p3_gc350_10us` | 0.3 | 350 (`Gc_const=350 --allow-unused`) | broad smooth halo, no localization: d > 0.9 to 2.7 mm and d > 0.5 to 5.5 mm from the wall; crack energy 43 J/m |
| `section_h0p3_gauss_eta0p3_30us` | 0.3 | Gaussian 12 % (`fracture_gaussian.i`), run to 30 us (pulse + first reflections + 5 us of quench) | faint ring at the wall (d <= 0.25) plus four faint diagonal streaks (d ~ 0.2-0.3, r = 14-29 mm) created between 14 and 22 us where the waves reflected by two adjacent free faces superpose; no node above d = 0.5; the quench (from 25 us) cut the vibration energy from 492 to 154 J/m by 30 us |
| `section_h0p3_eta1_weibull_10us` | 1.0 | Weibull, shape 10, mean 2721 (`MultiApps/fracture/input_files=fracture_weibull.i`) | 8-9 radial cracks with d = 1 cores, 2-3.5 mm long (d > 0.9 to 2.94 mm, d > 0.5 to 3.53 mm from the wall); crack energy 107 J/m; compare `section_h0p3_eta1_uniform_vs_weibull_10us.png` |

Note (2026-10-21, later): the inputs now also write `vel_x`, `vel_y`, `vel_mag` and
the toughness field `Gc_var` to the Exodus file and `v_max` to the CSV; the runs
above were made before that change and only hold `d` and `pulse_load_aux`.
The Weibull variant's scale was corrected afterwards so that the mean equals
`Gc_const` (the `eta1_weibull` run above had a mean 5 % low).

## Discovery job 13293853 (2026-10-21): 2D section, Gaussian scatter, eta = 1, 20-pulse train

Snapshots pulled while running (`discovery_snapshot_damage.png`, `discovery_snapshot_velocity.png`):
- 10 us: eight radial cracks of about 3 mm around the borehole (the scatter-seeded pattern).
- 18 us: spall. The compressive pulse reflects as tension at the four free faces and opens
  one long crack parallel to each face, 3-4 mm inside it, plus four diagonal cracks from the
  corners where two reflections superpose. This is the 32.3 mm "extent" in the CSV.
- 58 us (end of window 1): quench done, vibration energy 2.9 J/m of the 1153 J/m injected.
- 85 us (pulse 2): the radial cracks grow to 6-8 mm, the spall planes re-activate.
Conclusion: eta = 1 is far above the experiment (no face spall was observed there); use the
`fracture_extent_core` / `damaged_volume_core` measures (r <= 15 mm) to read the radial cracks.

## Face-spall lever study (2026-10-21 evening): eta = 0.5, Gaussian scatter, one window to 30 us

| run | lever | hole cracks (d >= 0.5, from the wall) | faces / corners after the reflections |
| --- | --- | --- | --- |
| `section_h0p3_gauss_eta0p5_baseline_14us.png` (cancelled at 14.8 us) | none | 7 cracks, 2.65 mm by 12 us | face-parallel spall bands starting at 14 us, 4-5 mm inside each face |
| `section_h0p3_gauss_eta0p5_attenuation_30us` | `zeta_pulse = 1e-8` s (Kelvin-Voigt loss factor ~0.06 at 1 MHz) | identical: 7 cracks, 2.65 mm, crack energy 47.5 J/m at 14 us | no face bands at all; only faint diagonal corner streaks (d <= 0.24) at 18-22 us; crack energy 58 J/m at 30 us |

| `section_h0p3_gauss_eta0p5_visc1_30us` | viscous damage, eta_v = 1 Pa s (`fracture_gaussian_viscous.i`) | smaller: 1.37 mm (r = 5.5 mm), crack energy 39 J/m at 14 us | corners still crack at 22 us (total extent 23.1 mm, 311 nodes d > 0.5 outside the core): worse on both counts |
| `section_h0p3_gauss_eta0p5_rise1us_30us` | rise time 1 us (shape 5.379e5 / 1.672e6), same 6 us duration | larger: 3.7 mm by 14 us, 10.8 mm by 22 us (the fuller pulse injects more impulse) | massive: faces and corners, total extent 32.3 mm, 1128 nodes d > 0.5 outside the core |
| `section_h0p3_gauss_eta0p5_zp1em7` (see CSV) | `zeta_pulse = 1e-7` s | same pattern as 1e-8 with ~6 % less crack energy | (run to 30 us; no benefit over 1e-8) |

Reading: the sharp front carries the spall tension; damping it at the material's own
1 MHz loss factor (zeta_pulse = 1e-8 s) removes the face spall without touching the radial
cracks. Rate-dependent damage slows the hole cracks but not the corner damage; a smoother
pulse at the same peak carries more impulse and makes everything worse. Picture:
`section_h0p3_gauss_eta0p5_levers_30us.png`.

## Discovery, 20-pulse production cases (eta = 0.5, zeta_pulse = 1e-8 s, Gaussian scatter), snapshot after two pulses

| case | hole cracks after pulse 1 / 2 (from the wall) | far field after pulse 2 |
| --- | --- | --- |
| free faces (job 13297781) | 2.65 / 6.48 mm, 7 cracks | four corner diagonals broken 20 us into pulse 2 (the pulse-1 streaks went critical); no face spall |
| absorbing faces (job 13297782) | 2.65 / 6.39 mm, 7 cracks | nothing (no node above d = 0.1); dashpots absorb 76 % of the pressure work |

Picture: `discovery_att_eta0p5_free_vs_absorb_pulse2.png`. Both jobs continue toward 20 pulses.
