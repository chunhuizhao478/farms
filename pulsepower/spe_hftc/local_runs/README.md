# Local 2D section runs (Mac, 8 ranks, JIT shim on PATH)

Section mesh `spe_section_h0p3` (l = 0.6 mm, Gc = 2721 J/m^2 from ft = 69 MPa),
one pulse, first 10 us only (no quench yet). Command:

    PATH=../tools/jit_shim:$PATH mpiexec -n 8 farms-opt -i elasticity.i Np=1 Executioner/end_time=1e-5 Outputs/file_base=<name> [eta=1.0 | Gc_const=350 --allow-unused | MultiApps/fracture/input_files=fracture_weibull.i | MultiApps/fracture/input_files=fracture_gaussian.i]

| run | eta | Gc (J/m^2) | result |
| --- | --- | --- | --- |
| `section_h0p3_eta0p3_10us` | 0.3 | 2721 | thin uniform damage ring, d_max = 0.18 within 0.55 mm of the wall; crack energy 5.6 J/m; no node above d = 0.5 |
| `section_h0p3_eta1_10us` | 1.0 | 2721 | crushed collar: d_max = 0.97, d > 0.5 out to 1.76 mm from the wall, d > 0.9 to 0.75 mm; crack energy 77 J/m; two incipient radial cracks (1-1.5 mm) break the symmetry |
| `section_h0p3_gc350_10us` | 0.3 | 350 (`Gc_const=350 --allow-unused`) | broad smooth halo, no localization: d > 0.9 to 2.7 mm and d > 0.5 to 5.5 mm from the wall; crack energy 43 J/m |
| `section_h0p3_eta1_weibull_10us` | 1.0 | Weibull, shape 10, mean 2721 (`MultiApps/fracture/input_files=fracture_weibull.i`) | 8-9 radial cracks with d = 1 cores, 2-3.5 mm long (d > 0.9 to 2.94 mm, d > 0.5 to 3.53 mm from the wall); crack energy 107 J/m; compare `section_h0p3_eta1_uniform_vs_weibull_10us.png` |

Note (2026-10-21, later): the inputs now also write `vel_x`, `vel_y`, `vel_mag` and
the toughness field `Gc_var` to the Exodus file and `v_max` to the CSV; the runs
above were made before that change and only hold `d` and `pulse_load_aux`.
The Weibull variant's scale was corrected afterwards so that the mean equals
`Gc_const` (the `eta1_weibull` run above had a mean 5 % low).
