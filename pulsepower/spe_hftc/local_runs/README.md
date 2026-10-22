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

## Attenuated 20-pulse cases, progress after 13 pulses (2026-10-22, 09:20 PDT)

Vista job 1089067 (72 ranks per case) and Discovery jobs 13297781/13297782 run the same two
inputs and agree to the step; the table is from the Vista CSVs. Crack extent from the borehole
wall at the end of each 60 us window, measured as the farthest node with d >= 0.5
(`fracture_extent`; the d >= 0.9 extent trails it by only 0.5 mm, so the cracks are fully
formed). `fracture_extent_core` saturates at 10.85 mm (15 mm core minus the 4.15 mm hole
radius) and is not used here. The free case loses its meaning after pulse 3: the corner
diagonals cut through, the near field crushes and the window maxima of v_max reach 1.1 km/s
by pulse 6 (fragments), so only its first three pulses are listed.

| after pulse | absorbing faces (mm) | growth (mm) | crack energy (J/m) | dashpot work (J/m) | free faces (mm) | paper (mm) |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | 2.65 | 2.65 | 47.5 | 537 | 2.65 | |
| 2 | 6.39 | 3.74 | 98.0 | 1078 | 6.48 (corner diagonals break) | |
| 3 | 8.35 | 1.96 | 127.4 | 1634 | 10.81 (corners cut through) | |
| 4 | 9.83 | 1.48 | 150.2 | 2191 | | |
| 5 | 11.04 | 1.21 | 168.4 | 2747 | | 7.97 |
| 6 | 12.42 | 1.38 | 183.2 | 3298 | | |
| 7 | 13.60 | 1.18 | 195.2 | 3844 | | |
| 8 | 14.49 | 0.89 | 205.3 | 4384 | | |
| 9 | 15.38 | 0.89 | 214.4 | 4919 | | |
| 10 | 16.27 | 0.89 | 222.9 | 5451 | | |
| 11 | 17.16 | 0.89 | 231.5 | 5978 | | 13.90 |
| 12 | 17.87 | 0.71 | 240.1 | 6501 | | |
| 13 | 18.65 | 0.78 | 248.9 | 7019 | | |
| 15 | | | | | | 15.66 |
| 20 | | | | | | 20.45 |

Reading. The pressure work is 695-734 J/m per pulse (rising slowly as the borehole gets more
compliant); the dashpots take 518-557 J/m of it per pulse, the cracks about 10 J/m. Seven
radial cracks grow all the way (picture `vista_att_eta0p5_damage_absorb788us_free386us.png`);
the longest one, toward +x, is about 3 mm from the side face after 13 pulses and reaches it
around pulse 17. Against the measurements (`vista_att_eta0p5_extent_vs_pulse.png`) the model
runs 23 % to 39 % long: 11.0 vs 7.97 mm at five pulses and 17.2 vs 13.9 mm at eleven. The
offset is made in the first five pulses (2.65, 3.74, 1.96, 1.48, 1.21 mm); from pulse 7 on the
model grows 0.7-0.9 mm per pulse, the same slope as the measured 0.83 mm per pulse between
pulses 5 and 20. The free case at 386 us is a fragmented sample: four corner diagonals with
branched tips, a crushed near field, 655 J/m of crack energy and 289 mm2 of damaged area.

Discovery jobs 13297781 (free) and 13297782 (absorbing) were cancelled on 2026-10-22 at 09:53 PDT
after 10 h 55 m, at 245 us (pulse 5) and 399 us (pulse 7, 13.60 mm, identical to Vista); the Vista
job 1089067 is the only run of these cases from here on.
The Vista free case was stopped by the user at 389 us (pulse 7) on 2026-10-22 at 10:00 PDT; the absorbing
case continues alone in job 1089067.

## Homogeneous Gc with eta = 0.4, absorbing faces: one pulse on the Mac (2026-10-22, 12:30 PDT)

Test asked for on 2026-10-22: the production design with the Gaussian Gc scatter removed
(`fracture.i`, constant Gc = 2721 J/m^2) and eta lowered from 0.5 to 0.4, absorbing side faces,
zeta_pulse = 1e-8 s, one 60 us window on 8 ranks (`section_h0p3_homog_eta0p4_absorb_1pulse.*`,
21 s per step). It is compared with the Vista Gaussian eta = 0.5 absorbing case at the same
instant in `section_h0p3_homog_eta0p4_absorb_vs_vista_eta0p5_10us.png` (40 mm view) and
`section_h0p3_homog_eta0p4_absorb_vs_vista_eta0p5_10us_zoom8mm.png` (16 mm close-up).

| at 10 us                               | homogeneous Gc, eta 0.4 (Mac) | Gaussian Gc, eta 0.5 (Vista) |
| -------------------------------------- | ----------------------------- | ---------------------------- |
| peak damage                            | 0.32, a uniform ring at the wall | 1.0                       |
| d >= 0.5 within 0.5 mm of the wall     | 0 %                           | 33 %                         |
| d >= 0.5 at 1-2 mm from the wall       | 0 %                           | 11 %                         |
| cracks longer than 1 mm                | none                          | 18, longest 2.4 mm           |
| pressure work                          | 594 J/m                       | 699 J/m                      |
| crack energy                           | 14.5 J/m (diffuse ring)       | 46 J/m                       |
| v_max                                  | 30 m/s                        | 41 m/s                       |

The homogeneous eta = 0.4 medium does not nucleate anything. The ring reaches d = 0.32 at 6 us
and is identical at 10 us: the pressure is below 10 % of its peak after 6 us and the absorbing
faces return nothing, so the remaining 50 us of the pulse cannot change it. Lowering eta from
0.5 to 0.4 cuts the peak wall pressure by 7 % (206 vs 222 MPa), and without weak spots the
driving force stays below the AT1 threshold everywhere, in line with the 2026-10-21 local tests
(eta 0.3 gave a faint ring, eta 1 a crushed collar, and only the scattered Gc localized radial
cracks). Taken together: the scatter is what nucleates the near-wall cracks, and eta = 0.4 sits
below the nucleation threshold of the homogeneous model. The two changes work against the
goal of more near-borehole fracture; a denser crushed zone needs a stronger near-wall driving
force (compressive split, smaller l, or borehole-fluid pressure in the cracks), not a weaker one.
