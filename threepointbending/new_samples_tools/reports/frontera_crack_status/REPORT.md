---
title: "New 3PB samples (c)-(f): crack status of the 16 Frontera restart runs"
subtitle: "Snapshot of 2026-10-15, 09:33 PDT; runs still going until their 24 h limits"
date: 2026-10-15
geometry: margin=1in
fontsize: 11pt
colorlinks: true
header-includes: |
  \usepackage{amsmath,amssymb,booktabs,siunitx}
---

## Summary

**Problem.** Sixteen new three-point-bending samples in four classes, (c) two holes, (d) connected holes,
(e) sine-wave void and (f) pillars, ran on Frontera to $t = 20$ and were restarted with a fine time step
($\Delta t = 0.001$) from 1.5 before each estimated crack time.
**Outcome.** All 16 cracked, at $t = 21.61$ to $29.73$, every one inside its fine window. The four cases also run on
USC Discovery cracked at the same time to within 0.002, so the coarse early stepping does not move the crack. At
this snapshot the crack tips sit at $y = 4.2$ to $7.4$ mm; none has reached the top (8.0 mm; pillars 8.5 mm), and D7
and D8 are held at the slot joining their holes.
**Risk.** The runs stop at their 24 h limits today, when Frontera's queues close.
**Out of scope.** Load-displacement curves and the experiments.

## Data Flow

1. Meshes, decks and jobs for the 16 STL cases, listed in the [case catalogue](#case-catalogue).
2. Production runs to $t = 20$ on Frontera (400 ranks): $\Delta t = 0.1$ to $t = 14$, then $0.001$.
3. Crack-time estimates: five $\Delta t = 0.1$ probe runs, a fit for the other eleven, see [Formulation](#formulation).
4. Restarts from the production checkpoints with the fine window \eqref{eq:window}.
5. Crack onsets read from the restart logs \eqref{eq:onset}, see [Results](#results).
6. Mid-thickness slices of the latest output, see [Crack snapshots](#crack-snapshots).

## Notation

| Symbol | Meaning | Code or file |
|--------|------------------------------|--------------------------|
| $t$ | load parameter; load-line displacement is 0.1 mm per unit $t$ | `func_loading = -1e-4 * t` |
| $d$ | phase-field damage, 0 intact, 1 broken | nodal field `d` |
| $\Delta t$ | time step | `dt_fn` in `elasticity.i` |
| $t_\mathrm{est}$ | estimated crack time, probe or fit | `logs/restart_schedule.tsv` |
| $t_\mathrm{fine}$ | start of the $\Delta t = 0.001$ window | `make_restart_decks.py` |
| $N_\mathrm{fp}(t_n)$ | fixed-point iterations of step $n$ | restart log, `Beginning fixed point iteration` |
| $t_\mathrm{on}$ | crack onset, first hard step | `scripts/frontera_hard_steps.sh` |
| $y_\mathrm{tip}$ | crack-tip height, whole mesh | `scripts/frontera_extract_slices.py` |
| $y_\mathrm{top}$ | beam top, the load line: 8.0 mm, pillars 8.5 mm | mesh |

## Formulation

\begin{equation}\label{eq:fit}
t_\mathrm{est} = 14.40 + \frac{1.521}{\max d\,(t = 15.5)}
\end{equation}
*Fitted* - `FIT_A`, `FIT_B` in `make_restart_decks.py`, from the five probe crack times, see [Appendix A.1](#a.1).

- $t_\mathrm{est}$: estimated crack time (load units).
- $\max d\,(t = 15.5)$: largest damage in the production run at $t = 15.5$ (dimensionless).

\begin{equation}\label{eq:window}
t_\mathrm{fine} = \max\Bigl(20,\ \tfrac{1}{2}\bigl\lfloor 2\,(t_\mathrm{est} - 1.5)\bigr\rfloor\Bigr), \qquad
\Delta t = \begin{cases} 0.1 & t < t_\mathrm{fine} \\ 0.001 & t \ge t_\mathrm{fine} \end{cases}
\end{equation}
*Definition* - `make_restart_decks.py`: lead 1.5, rounded down to a multiple of 0.5.

- $t_\mathrm{fine}$: start of the fine window (load units).
- $\Delta t$: time step of the restart (load units).

\begin{equation}\label{eq:onset}
t_\mathrm{on} = \min\{\, t_n : N_\mathrm{fp}(t_n) \ge 5 \,\}
\end{equation}
*Definition* - the staggered elasticity-damage loop needs many iterations when a crack jumps; one iteration otherwise.

- $t_n$: time of step $n$; $N_\mathrm{fp}(t_n)$: its fixed-point iterations (cap 20).

\begin{equation}\label{eq:tip}
y_\mathrm{tip} = \max\{\, y_i : d_i > 0.95 \,\}
\end{equation}
*Definition* - over all mesh nodes $i$ at the last complete output.

- $y_i$: node height (mm); $d_i$: its damage.

## Called Code

| Script | Role | Runs on |
|----------------------------|------------------------------|------------|
| `make_restart_decks.py` | restart decks with the fine window | Mac |
| `scripts/frontera_hard_steps.sh` | fixed-point iterations per step, crack onsets | Frontera, read-only |
| `scripts/frontera_extract_slices.py` | tetrahedra crossing $z = 4$ mm, latest output | Frontera, read-only |
| `render_figures.py` | exact cut at $z = 4$ mm, figures, tables | Mac |
| `build.sh` | this PDF | Mac |

## Phases at a Glance

| Phase | What | Status |
|---|---|---|
| 1 | meshes, decks and jobs for 16 cases | done 2026-10-10 |
| 2 | production runs to $t = 20$ | done 2026-10-13 |
| 3 | five $\Delta t = 0.1$ probe runs | done 2026-10-12 |
| 4 | restarts with the fine window | running; limits end 2026-10-15 |
| 5 | cross-check: four cases on USC Discovery | done 2026-10-14 |

## Case Catalogue {#case-catalogue}

IDs used in every table and panel. Lengths in mm; $x$ along the beam (notch at $x = 14$), $y$ up (notch tip at
$y = 1.6$), $z$ through the 8 mm thickness.

| ID | Class | Case | Geometry (mm) | STL file | Frontera job |
|---|----|----------------------|--------------------------|------------------------|--------|
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

## Results {#results}

Onset from \eqref{eq:onset}; estimate from a probe run or \eqref{eq:fit}, and the same case on USC Discovery where
it ran; tip from \eqref{eq:tip} at the snapshot output; "to top" is $y_\mathrm{top} - y_\mathrm{tip}$.

| ID | Onset t | Estimate t (source) | Fine dt from | Snapshot t | Tip y (mm) | To top (mm) | Tip position |
|---|-------|------------------|------|-------|------|------|----------------------|
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

- Every onset falls inside its fine window; onsets come from 0.1 after to 0.7 before the estimate.
- The four cases run on both machines (C1, D2, E1, F1) agree to within 0.002 in $t_\mathrm{on}$.
- D7 and D8 are held at the bottom of their slot ($y = 4.20$); D3 and D9 sit just below theirs.
- D2 on Discovery waited at its slot until that run was cancelled at $t = 28.7$; here it passed and reached 7.32 mm.
- In (e) the crack runs through the void where it crosses the crack line (E1) and passes it otherwise (E2).
- In (f) the crack crosses the pillar gap and continues about 0.5 mm to the left of the lower crack.
- Growth is slow over the last millimetre, under the load line; no crack has reached the top yet.

## Crack Snapshots {#crack-snapshots}

Each panel shows the damage $d$ on the plane $z = 4$ mm (mid-thickness) at the snapshot output, cut exactly through
the tetrahedra. White is void (notch, holes, slot, sine void, gaps between pillars) or outside the beam. The cyan
triangle marks the highest broken point ($d > 0.95$) in this plane; the label gives the whole-mesh tip
\eqref{eq:tip}. The window is $x = 7$ to $21$ mm around the notch.

![(c) two holes, not connected: C1 to C3.](figures/fig_c_asym.png){width=100%}

![(d) connected holes, left hole D = 1.5 mm: D1 to D5.](figures/fig_d_conn_DL15.png){width=100%}

![(d) connected holes, both D = 2.0 mm: D6 to D9.](figures/fig_d_conn_DL20.png){width=100%}

![(e) sine-wave void: E1 centred, E2 offset by one amplitude.](figures/fig_e_sine.png){width=100%}

![(f) pillar structures, welded mismatch: F1 2onN, F2 3onN.](figures/fig_f_pillars.png){width=100%}

## Sources

- [S1] STL files: `Paper_ThreePointBending/new_samples/New sample stls/`, one folder per class.
- [S2] Frontera outputs: `$SCRATCH2/projects/farms_cdms_04192026/threepointbending/code_<family>/case_*/`.
- [S3] Frontera restart logs: `jobs_<family>_recover/3dptrbend_*_recover.o<job>`; job map in `data/job_case_map.txt`.
- [S4] Crack onsets: `data/hard_steps.txt`, scanned 2026-10-15 09:37 PDT.
- [S5] USC Discovery runs: condo jobs 6406191 (C1), 6406822 (D2), 6413383 (E1), 6439163 (F1).
- [S6] Code: farms_cdms branch `worktree-tpb-new-samples`, `threepointbending/new_samples_tools/`.

## Appendix A: Derivations

### A.1 Crack-time fit \eqref{eq:fit} {#a.1}

The five probe runs ($\Delta t = 0.1$ to $t = 50$) gave crack times 21.8, 22.9, 24.5, 28.0 and 30.0 (C1, E1, F1,
D2, D9). Against the largest damage of the same cases' production runs at $t = 15.5$ (0.2011, 0.1760, 0.1550, 0.1136,
0.0958) they follow $t_\mathrm{est} = a + b / \max d$. The decks use $a = 14.40$, $b = 1.521$, with residuals of at
most 0.3 on the five points; an ordinary least-squares fit gives $a = 14.45$, $b = 1.512$, which moves every
$t_\mathrm{est}$ by less than 0.05. The restart onsets fall from 0.1 after to 0.7 before the estimates, inside the
1.5 lead of \eqref{eq:window}.

## Appendix B: Rationale

None yet.
