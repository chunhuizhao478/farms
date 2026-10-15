| ID | Class | Case | Geometry (mm) | STL file | Frontera job |
|---|---|---|---|---|---|
| C1 | (c) | `xL3.00mm_xR2.00mm_y3.6mm_radi1.0mm` | holes r 1.0 at x 14-3.00 and 14+2.00, y 5.2 | `2mm_circle_3.6up_2.0-3.0mm apart v2.stl` | 7990425 |
| C2 | (c) | `xL3.00mm_xR2.50mm_y3.6mm_radi1.0mm` | holes r 1.0 at x 14-3.00 and 14+2.50, y 5.2 | `2mm_circle_3.6up_2.5-3.0mm apart v2.stl` | 7990426 |
| C3 | (c) | `xL3.50mm_xR2.50mm_y3.6mm_radi1.0mm` | holes r 1.0 at x 14-3.50 and 14+2.50, y 5.2 | `2mm_circle_3.6up_2.5-3.5mm apart v2.stl` | 7990427 |
| D1 | (d) | `DL1.5mm_DR2.0mm_x1.00mm_y3.6mm_overlap` | holes r 0.75/1.00 at x 14 -/+ 1.00, y 5.2; holes touch, no slot | `Alum_2.0_1.5mmD_y3.6mm_x1.0mm_overlap.stl` | 7990428 |
| D2 | (d) | `DL1.5mm_DR2.0mm_x1.50mm_y3.6mm_downline` | holes r 0.75/1.00 at x 14 -/+ 1.50, y 5.2; slot y 4.45-4.65 | `Alum_2.0_1.5mmD_y3.6mm_x1.5mm_downline.stl` | 7990429 |
| D3 | (d) | `DL1.5mm_DR2.0mm_x1.50mm_y3.6mm_upline` | holes r 0.75/1.00 at x 14 -/+ 1.50, y 5.2; slot y 5.75-5.95 | `Alum_2.0_1.5mmD_y3.6mm_x1.5mm_upline.stl` | 7990430 |
| D4 | (d) | `DL1.5mm_DR2.0mm_x2.00mm_y3.6mm_downline` | holes r 0.75/1.00 at x 14 -/+ 2.00, y 5.2; slot y 4.45-4.65 | `Alum_2.0_1.5mmD_y3.6mm_x2.0mm_downline.stl` | 7990431 |
| D5 | (d) | `DL1.5mm_DR2.0mm_x2.00mm_y3.6mm_midline` | holes r 0.75/1.00 at x 14 -/+ 2.00, y 5.2; slot y 5.10-5.30 | `Alum_2.0_1.5mmD_y3.6mm_x2.0mm_midline.stl` | 7990432 |
| D6 | (d) | `DL2.0mm_DR2.0mm_x1.00mm_y3.6mm_overlap` | holes r 1.00/1.00 at x 14 -/+ 1.00, y 5.2; holes touch, no slot | `Alum_2.0_2.0mmD_y3.6mm_x1.0mm_overlap.stl` | 7990433 |
| D7 | (d) | `DL2.0mm_DR2.0mm_x1.50mm_y3.6mm_downline` | holes r 1.00/1.00 at x 14 -/+ 1.50, y 5.2; slot y 4.20-4.40 | `Alum_2.0_2.0mmD_y3.6mm_x1.5mm_downline.stl` | 7990434 |
| D8 | (d) | `DL2.0mm_DR2.0mm_x2.00mm_y3.6mm_downline` | holes r 1.00/1.00 at x 14 -/+ 2.00, y 5.2; slot y 4.20-4.40 | `Alum_2.0_2.0mmD_y3.6mm_x2.0mm_downline.stl` | 7990435 |
| D9 | (d) | `DL2.0mm_DR2.0mm_x2.00mm_y3.6mm_upline` | holes r 1.00/1.00 at x 14 -/+ 2.00, y 5.2; slot y 6.00-6.20 | `Alum_2.0_2.0mmD_y3.6mm_x2.0mm_upline.stl` | 7990436 |
| E1 | (e) | `D1.0mm_y1.6mm_amp2.0mm_centered` | sine tube r 0.5 at y 3.2, amplitude 2, axis x 14.0 | `Alum_1mmD_1.6up_sinewavevoid_Amp2_centered.stl` | 7990438 |
| E2 | (e) | `D1.0mm_y1.6mm_amp2.0mm_1ampoffset` | sine tube r 0.5 at y 3.2, amplitude 2, axis x 16.0 | `Alum_1mmD_1.6up_sinewavevoid_Amp2_1ampoffset.stl` | 7990437 |
| F1 | (f) | `welded_mismatch_2onN` | blocks gap y 4.0-4.5, pillars r 0.5, 2onN rows | `Welded_mismatch_2onN v2.stl` | 7990439 |
| F2 | (f) | `welded_mismatch_3onN` | blocks gap y 4.0-4.5, pillars r 0.5, 3onN rows | `Welded_mismatch_3onN v2.stl` | 7990440 |

| ID | Onset t | Estimate t (source) | Fine dt from | Frame t | Tip y (mm) | To top (mm) | Tip position |
|---|---|---|---|---|---|---|---|
| C1 | 21.614 | 21.8 (probe); Discovery 21.614 | 20.0 | 24.7 | 6.42 | 1.58 | above the hole tops (y 6.2) |
| C2 | 22.392 | 22.4 (fit) | 20.5 | 25.1 | 6.55 | 1.45 | above the hole tops (y 6.2) |
| C3 | 23.191 | 23.1 (fit) | 21.5 | 23.6 | 6.26 | 1.74 | above the hole tops (y 6.2) |
| D1 | 24.308 | 24.5 (fit) | 23.0 | 32.2 | 7.27 | 0.73 | above the touching holes |
| D2 | 27.749 | 28.0 (probe); Discovery 27.749 | 26.5 | 33.8 | 7.32 | 0.68 | past the slot |
| D3 | 26.966 | 27.5 (fit) | 25.5 | 32.6 | 5.60 | 2.40 | below the slot |
| D4 | 27.897 | 28.2 (fit) | 26.5 | 32.2 | 7.15 | 0.85 | past the slot |
| D5 | 28.796 | 29.5 (fit) | 27.5 | 32.9 | 7.13 | 0.87 | past the slot |
| D6 | 27.954 | 28.1 (fit) | 26.5 | 35.4 | 7.31 | 0.69 | above the touching holes |
| D7 | 25.876 | 25.8 (fit) | 24.0 | 31.4 | 4.20 | 3.80 | held at the slot bottom |
| D8 | 27.941 | 28.3 (fit) | 26.5 | 31.9 | 4.20 | 3.80 | held at the slot bottom |
| D9 | 29.731 | 30.0 (probe) | 28.5 | 33.6 | 5.82 | 2.18 | below the slot |
| E1 | 22.719 | 22.9 (probe); Discovery 22.720 | 21.0 | 30.1 | 7.01 | 0.99 | past the void band (y 2.7-3.7) |
| E2 | 22.692 | 22.9 (fit) | 21.0 | 29.8 | 7.07 | 0.93 | past the void band (y 2.7-3.7) |
| F1 | 24.242 | 24.5 (probe); Discovery 24.242 | 23.0 | 27.3 | 7.37 | 1.13 | past the pillar gap (y 4.0-4.5) |
| F2 | 24.506 | 24.4 (fit) | 22.5 | 26.8 | 7.25 | 1.25 | past the pillar gap (y 4.0-4.5) |
