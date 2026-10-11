---
title: "New 3PB samples (c)-(f): meshes, decks and jobs"
date: 2026-10-10
geometry: margin=1in
fontsize: 11pt
colorlinks: true
header-includes: |
  \usepackage{amsmath,amssymb,booktabs,siunitx}
---

## Summary

**The problem.** Sixteen new specimens exist only as STL drawings, in four families: two holes at unequal offsets, hole pairs joined by a thin line, a sine-wave void, and blocks joined by pillars. The phase-field runs need meshes, decks and jobs like the earlier sets.

**The fix.** Each STL is measured, rebuilt as exact CAD geometry in the two-hole mesh layout, meshed behind a memory kill switch, and checked against its STL. The two-hole deck and job are reused with only the mesh path changed.

**Expected outcome.** 16 meshes of 1.3-2.1 million tetrahedra at the two-hole resolution, each within 5 um of its STL (sine voids: 17 um). MOOSE accepts all 16 decks.

**Main risk / tradeoff.** The refined size, 100 um, is twice the phase-field length. The notch tip stays round; the STLs draw it flat.

**Out of scope.** Running the jobs: Frontera's farms build was purged.

## Data Flow

Nine stages, in the order the code runs. A bare name lives in memory, a path is a file, and **bold** marks what ships. Paths are relative to `threepointbending/`; `<family>` is `2holes_asym`, `2holes_connected`, `sinewave` or `pillars`.

**Geometry**

1. **Read the STLs.** Sixteen binary STLs in mm; circle and sine fits give each case's numbers, kept as `CASES`. [Phase\ 1](#phase-1)
2. **Write the geometry.** One OpenCASCADE script per case, built from \eqref{eq:holes}-\eqref{eq:sine}: `meshfile_<family>/<mesh>.geo`. [Phase\ 1](#phase-1)
3. **Check the geometry.** Gmsh builds each solid, whose volume must equal \eqref{eq:volume}: report. [Phase\ 1](#phase-1)

**Mesh**

4. **Mesh it.** Gmsh, one case at a time behind the memory kill switch: **`meshfile_<family>/<mesh>.msh`** and `new_samples_tools/logs/mesh_log.tsv`. [Phase\ 2](#phase-2)
5. **Compare with the STL.** Two-way surface distance, volume, line-support nodes and quality \eqref{eq:quality}: `new_samples_tools/logs/verify_log.tsv`. [Phase\ 3](#phase-3)
6. **Export for ParaView.** A legacy VTK copy beside each mesh: **`meshfile_<family>/<mesh>.vtk`**. [Phase\ 5](#phase-5)

**Simulation inputs**

7. **Write decks and jobs.** The two-hole deck and job with the mesh path changed: **`code_<family>/case_<family>_3d_<tag>/`**, **`jobs_<family>/`**, **`jobs_<family>_recover/`**. [Phase\ 4](#phase-4)
8. **Check the decks with MOOSE.** Every deck through `--check-input`, the largest through the jobs' serial `--split-mesh 400`: `new_samples_tools/logs/check_decks_moose.tsv`. [Phase\ 5](#phase-5)

**Runs (not started)**

9. **Stage and submit on Frontera.** Restore the purged farms tree, rebuild `farms-opt`, upload the meshes, submit: job outputs. [Phase\ 6](#phase-6)

Stages 2-5 build exact geometry and prove it against the STLs; stage 8 proves each deck reads its mesh with the MOOSE commit that Frontera uses.

## Notation

| Symbol | Meaning | Units | Code variable |
| ------ | ------- | ----- | ------------- |
| $x, y, z$ | length, height and thickness coordinates | m | `x`, `y`, `z` |
| $L, H, T$ | specimen length 28, height 8, thickness 8 mm | m | `BEAM_*_MM`, `extrude_z` |
| $x_N, y_N$ | notch tip, (14, 1.6) mm | m | `notch_tip_x`, `notch_tip_y` |
| $w_N, d_N$ | notch width 1 mm and depth 1.6 mm | mm | `NOTCH_WIDTH_MM`, `NOTCH_TIP_MM` |
| $A_N$ | notch cross-section area | mm$^2$ | `notch_area_mm2` |
| $\Delta x_L, \Delta x_R$ | hole offsets left and right of the notch tip | m | `dx_left`, `dx_right` |
| $\Delta y$ | hole or void height above the notch tip | m | `dy` |
| $x_L, x_R, y_c$ | hole centres | m | `hole_left_x`, `hole_right_x`, `hole_y` |
| $r, r_L, r_R$ | hole radii; $r_{\min}$ is the smaller | m | `hole_radius_left`, `hole_radius_right` |
| $w$ | connecting-line width, 0.2 mm | m | `line_width` |
| $y_0, y_1$ | lower and upper edge of a line or throat | mm | `connected_band_mm` |
| $c_L, c_R$ | half chords of the holes at height $y$ | mm | `half_chord` |
| $w_c$ | narrowest material cusp kept, 0.1 mm | m | `cusp_min_width` |
| $s_0$ | throat half height | m | `throat_h` |
| $x_0, A, \lambda$ | sine axis, amplitude 2 mm, period 8 mm | m | `x_axis`, `amp`, `wavelength` |
| $r_v, y_v$ | sine-void radius 0.5 mm and centre height | m | `void_radius`, `void_y` |
| $n_p, r_p$ | pillar count and radius 0.5 mm | -, m | `row_n`, `pillar_radius` |
| $V, V_{\mathrm{void}}$ | solid volume, void volume | mm$^3$ | `expected_volume_mm3` |
| $V_e, \ell_k$ | tetrahedron volume and its six edge lengths | m$^3$, m | `vol`, `edges` |
| $\eta$ | mean-ratio quality, 1 for a regular tetrahedron | - | `q` |
| $h_r, h_s, h_g$ | refined, support and global element size | m | `lc_refined`, `lc_support`, `lc` |
| $l$ | phase-field length, 50 um | m | `l` in `elasticity.i` |
| $t, dt$ | deck time and time step (0.1, then 0.001 after $t$ = 14) | - | `t`, `dt_fn` |

## Formulation

\begin{equation}\label{eq:holes}
x_L = x_N - \Delta x_L, \qquad x_R = x_N + \Delta x_R, \qquad y_c = y_N + \Delta y
\end{equation}
*Definition* - offsets measured from the STLs [S1](#sources); a wrong offset moves a hole.

\begin{equation}\label{eq:line}
[y_0, y_1] = \begin{cases}
[y_c + r_{\min} - w,\; y_c + r_{\min}] & \text{up line} \\
[y_c - w/2,\; y_c + w/2] & \text{mid line} \\
[y_c - r_{\min},\; y_c - r_{\min} + w] & \text{down line}
\end{cases}
\end{equation}
*Definition* - measured from the STLs [S1](#sources): up and down lines are tangent to the smaller hole.

\begin{equation}\label{eq:throat}
s_0 = \sqrt{r^2 - (r - w_c/2)^2}, \qquad [y_0, y_1] = [y_c - s_0,\; y_c + s_0]
\end{equation}
*Derived* - see [Appendix A.1](#a.1). It joins touching holes, $x_R - x_L = 2r$, without unmeshable cusps.

\begin{equation}\label{eq:sine}
\bigl(x - x_0 - A \sin(2 \pi z / \lambda)\bigr)^2 + (y - y_v)^2 \le r_v^2, \qquad y_v = y_N + \Delta y
\end{equation}
*Definition* - the void of family (e), fitted to the STL centreline within 0.03 mm [S12](#sources).

\begin{equation}\label{eq:volume}
V = T\,(L H - A_N) - V_{\mathrm{void}}, \qquad A_N = w_N (d_N - w_N/2) + \pi w_N^2 / 8
\end{equation}
*Derived* - see [Appendix A.2](#a.2), which also gives $V_{\mathrm{void}}$ per family.

\begin{equation}\label{eq:quality}
\eta = \frac{12\,(3 |V_e|)^{2/3}}{\sum_{k=1}^{6} \ell_k^2}
\end{equation}
*Definition* - the mean-ratio shape measure: 1 for a regular tetrahedron, 0 for a flat one.

- Element sizes: $h_r$ = 100 um in the feature boxes, $h_s$ = 250 um at the supports, $h_g$ = 1 mm elsewhere; the pillar layer gets 125 um.

## Approach

`make_new_sample_geos.py` holds `CASES` and writes each `.geo` in the `meshfile_2holes` layout: outline points 25/26/27 become the line-support edges, and Box fields set the sizes. Family (d) cuts its holes as one explicit outline, because tangent slot faces break OpenCASCADE booleans. Family (f) turns off union face-merging so the support edges survive. `mesh_new_samples.py` runs gmsh behind the kill switch. `make_new_sample_cases.py` copies the two-hole deck and jobs through asserted one-line substitutions.

## Called Code

| Function | What it computes | Governing equation | Verified how |
| -------- | ---------------- | ------------------ | ------------ |
| gmsh OCC `Cylinder`, `BooleanDifference` | beam minus two holes, family (c) | \eqref{eq:holes} | ran: volume equals \eqref{eq:volume} [S4](#sources) |
| gmsh OCC `ThruSections` | sine void as a loft through 49 circles | \eqref{eq:sine} | ran: volume within 3e-7 [S4](#sources) |
| gmsh OCC `BooleanUnion` | blocks plus pillars, family (f) | \eqref{eq:volume} | ran: volume and line nodes [S6](#sources) |
| gmsh `Box` and `Min` fields | element size targets | none - sizing | band median edge 131 um [S8](#sources) |
| `meshio.read` | MSH 2.2 tetrahedra and physical tags | none - I/O | ran on all 16 [S6](#sources) |
| VTK `compute_implicit_distance` | exact point-to-triangle distance | none - check | ran on all 16 [S6](#sources) |
| MOOSE mesh and nodeset generators | mesh, elastic pads, line nodesets | none - setup | `--check-input` on all 16 [S7](#sources) |
| `farms-opt --split-mesh 400` | serial partition used by every job | none - setup | ran on the largest mesh [S7](#sources) |

## Constraints

- Decks equal `code_2holes` except the mesh path; pillars also move the loading height to 8.5 mm.
- Meshes are MSH 2.2 with one Physical Volume, which MOOSE reads as block 1.
- Line supports need nodes exactly on $x$ = 4, 24 mm at the bottom and $x$ = 14 mm on top.
- Local runs stay below 20 GB RSS and above 12 percent free memory (36 GB laptop).
- Frontera queues close on 2026-10-15 [S10](#sources), and its farms tree is purged [S9](#sources).

## Phases at a Glance

| # | After this phase | Stages | Files touched | Done when |
| - | ---------------- | ------ | ------------- | --------- |
| [1](#phase-1) | every case has an exact, checked geometry | 1-3 | `make_new_sample_geos.py`, 16 `.geo` | 16 of 16 volumes match |
| [2](#phase-2) | every case has a mesh | 4 | `mesh_new_samples.py`, 16 `.msh` | the guard never fires |
| [3](#phase-3) | every mesh is proven against its STL | 5 | `verify_new_sample_meshes.py` | 16 of 16 PASS |
| [4](#phase-4) | every case has a deck and two jobs | 7 | `make_new_sample_cases.py`, `code_*`, `jobs_*` | diffs show only intended lines |
| [5](#phase-5) | MOOSE accepts every deck; ParaView files exist | 6, 8 | `check_decks_moose.py`, `export_vtk_for_paraview.py` | 16 of 16 checks pass |
| [6](#phase-6) | the 16 runs are queued on Frontera | 9 | Frontera project tree | jobs accepted by `sbatch` |

## Open Questions

| Question | What it blocks | How to resolve |
| -------- | -------------- | -------------- |
| Keep the round notch tip, or the STL's flat tip? | [Phase 2](#phase-2) meshes | user; set `NOTCH_ROUND = 0` and rerun |
| Keep 100 um refinement, or 50 um as the angles/height sets? | [Phase 2](#phase-2), run cost | user; about 8x elements |
| Run on Frontera before 2026-10-15, or on Discovery? | [Phase 6](#phase-6) | user |
| Is $dt$ = 0.1 up to $t$ = 14 still before crack onset? | [Phase 6](#phase-6) results | read the first run's damage history |
| Is the throat trim right for the touching holes? | one case | user check in ParaView |

## Phase 1: Geometry {#phase-1}

**In one sentence.** Turns each STL's measured numbers into an exact OpenCASCADE solid, written as one `.geo` per case.

**Implements [Data Flow](#data-flow) stages:** 1-3

**Files.** `new_samples_tools/make_new_sample_geos.py` (new): `CASES`, one writer per family, `--check`. `meshfile_<family>/mesh_<family>_3d_<tag>.geo` (new, 16 files).

**Requirements.**

1. Base beam as `meshfile_2holes`: outline points 1-9, 14, 25, 26, 27 extruded by $T$; switch `notch_round`.
2. Family (c): two `Cylinder`s placed by \eqref{eq:holes}, then `BooleanDifference`.
3. Family (d): one inner outline of the section with lines at $y_0, y_1$ from \eqref{eq:line} or \eqref{eq:throat}; arcs split at multiples of 90 degrees; no boolean.
4. Family (e): loft through 49 circles of \eqref{eq:sine} for $z$ from -0.5 to 8.5 mm, then `BooleanDifference`.
5. Family (f): notched lower block, upper block and $n_p$ pillars embedded 0.1 mm into both; `BooleanUnion` with `Geometry.OCCUnionUnify = 0`.
6. Fields 1/3/4/5/7 as `meshfile_2holes`; (e) one void box; (f) a gap box at $h_r$ and the full layer at 125 um.

**Edge cases.** A tangent line has a zero half chord, so the arc ends at the tangent point. Touching holes get the throat \eqref{eq:throat}.

**Acceptance criteria.**

- [x] Every `.geo` builds exactly one volume.
- [x] The OpenCASCADE volume equals \eqref{eq:volume} within relative 2e-5.

**Depends on:** nothing

## Phase 2: Mesh generation {#phase-2}

**In one sentence.** Runs gmsh on each `.geo`, one at a time, and kills it if memory runs short.

**Implements [Data Flow](#data-flow) stages:** 4

**Files.** `new_samples_tools/mesh_new_samples.py` (new). `meshfile_<family>/*.msh` (new, gitignored).

**Requirements.**

1. `run_guarded(cmd, log_path, cwd, env)` samples RSS (`ps`) and system free memory (`memory_pressure -Q`) every 2 s; SIGKILL above 20 GB or below 12 percent free.
2. `gmsh -3 <geo> -o <msh>`; MSH 2.2 is set inside the `.geo`.
3. Append nodes, tetrahedra, time, peak RSS and lowest free memory to `logs/mesh_log.tsv`.

**Acceptance criteria.**

- [x] 16 meshes written without the guard firing.
- [x] Each `.msh` holds only tetrahedra of one Physical Volume.

**Depends on:** [Phase 1](#phase-1)

## Phase 3: Verification against the STLs {#phase-3}

**In one sentence.** Reads each mesh with its STL and proves they describe the same solid, with nodes where the deck's line supports need them.

**Implements [Data Flow](#data-flow) stages:** 5

**Files.** `new_samples_tools/verify_new_sample_meshes.py` (new); `logs/verify_log.tsv`.

**Requirements.**

1. STL to model frame: mm to m; the sine STLs are rotated $(x, y, z) \to (x, z, 8\,\mathrm{mm} - y)$.
2. Distances both ways with VTK, excluding the two designed departures: the notch tip and the throat.
3. Each line support has at least 9 nodes spanning $z$ = 0 to 8 mm, and MOOSE's +-10 um boxes select exactly these.
4. Quality \eqref{eq:quality}, orientation, orphan nodes, and volume against \eqref{eq:volume}.

**Acceptance criteria.**

- [x] Both distances below 10 um, and below 40 um for the sine voids (STL spline versus sine).
- [x] $\eta > 0.05$ everywhere; no inverted tetrahedra; no orphan nodes; volume within relative 1e-3.

**Depends on:** [Phase 2](#phase-2)

## Phase 4: Decks and jobs {#phase-4}

**In one sentence.** Copies the two-hole deck and the Frontera jobs for each case, changing only names, the mesh path and, for pillars, the loading height.

**Implements [Data Flow](#data-flow) stages:** 7

**Files.** `new_samples_tools/make_new_sample_cases.py` (new). Per case: `code_<family>/case_<family>_3d_<tag>/elasticity.i` and `fracture.i`, `jobs_<family>/job_<tag>.sbatch`, `jobs_<family>_recover/job_<tag>_recover.sbatch`; one `submit_all.sh` per jobs folder.

**Requirements.**

1. Templates: deck `code_2holes/case_2holes_3d_x1.5mm_y1.6mm_radi0.5mm`, jobs `jobs_2holes` and `jobs_angles_recover`.
2. Every substitution matches exactly once (asserted).
3. Pillars: the elastic pad, loading-line box and centre node move from $y$ = 8 mm to 8.5 mm.
4. `submit_all.sh` changes into its own folder first, so `../../farms-opt` resolves.

**Acceptance criteria.**

- [x] `diff` against the templates shows only the mesh path, names and the four pillar loading-height lines.

**Depends on:** [Phase 1](#phase-1)

## Phase 5: MOOSE checks and ParaView export {#phase-5}

**In one sentence.** Feeds every deck to the local `farms-opt` without solving, measures the jobs' serial split step, and writes a `.vtk` per mesh.

**Implements [Data Flow](#data-flow) stages:** 6, 8

**Files.** `new_samples_tools/check_decks_moose.py` and `export_vtk_for_paraview.py` (new); `meshfile_<family>/*.vtk` (gitignored).

**Requirements.**

1. The local `farms-opt` is MOOSE `94b6528`, the commit of Frontera's `moose94b6528` build [S11](#sources).
2. `--check-input` per deck and `--split-mesh 400` on the largest mesh, both under `run_guarded`, in a working folder outside the repository.
3. `gmsh <msh> -save -format vtk -bin -o <vtk>`.

**Acceptance criteria.**

- [x] 16 of 16 `--check-input` runs exit 0 without `*** ERROR`.
- [x] The split's peak RSS stays below 10 percent of a Frontera node's 192 GB.
- [x] Every `.vtk` loads in VTK and holds only tetrahedra.

**Depends on:** [Phase 3](#phase-3), [Phase 4](#phase-4)

## Phase 6: Frontera staging and submission (pending) {#phase-6}

**In one sentence.** Restores the purged farms tree on Frontera, rebuilds `farms-opt`, uploads the meshes and queues the 16 jobs; nothing here has run.

**Implements [Data Flow](#data-flow) stages:** 9

**Files.** `new_samples_tools/frontera/build_farms_frontera.sbatch` (new, not submitted). Frontera `/scratch2/10024/zhaochun/projects/farms_cdms_04192026`, whose files were all purged [S9](#sources).

**Requirements.**

1. Check this branch out into that tree.
2. Build `farms-opt` with `build_farms_frontera.sbatch`: `gcc/9.1.0`, `mvapich2-x/2.3`, `MOOSE_DIR=$WORK/projects/moose94b6528` [S11](#sources).
3. Upload the 16 `.msh` files into `threepointbending/meshfile_<family>/`.
4. Run `submit_all.sh` in each `jobs_<family>` folder.

**Acceptance criteria.**

- [ ] `farms-opt --version` prints `94b6528` on a compute node.
- [ ] All 16 jobs are accepted by `sbatch`.

**Depends on:** [Phase 5](#phase-5) and the user's approval of every Frontera action.

## Testing Strategy

Geometry: the OpenCASCADE volume must equal the analytic volume. Meshes: two-way distance to the STL, volume, orientation, orphan nodes, the exact line-support nodes MOOSE selects, and mean-ratio quality. Decks: MOOSE `--check-input` on all 16 with the commit Frontera uses. Resolution: the refined band's median edge must equal that of the existing two-hole mesh built with the same setting.

## Risks

| Risk | Severity | Mitigation |
| ---- | -------- | ---------- |
| 100 um is twice the phase-field length; damage bands may smear | med | 50 um is one setting; about 8x elements |
| The round notch tip departs from the STL's flat tip | low | `NOTCH_ROUND = 0`, regenerate, rerun the checks |
| Frontera closes before the 48 h runs finish | high | submit by 2026-10-13, or use Discovery |
| The dt switch at $t$ = 14 may come after crack onset | med | inspect the first run's damage history |
| The throat changes the touching-holes crack path | low | it trims 0.02 mm$^2$; the X-ray shows merged holes |

## Sources

| Id | What it is | Where it lives |
| -- | ---------- | -------------- |
| S1 | 16 STLs in mm; circle fits leave residuals below 1e-6 mm | `Paper_ThreePointBending/new_samples/New sample stls` |
| S2 | Slide deck: families and experiment names | `new_samples/New Samples.pptx` |
| S3 | X-ray frames: round notch tip, merged figure-8 holes | `New Samples.pptx`, slide 10 |
| S4 | OpenCASCADE volumes equal analytic volumes, 16 of 16 | `make_new_sample_geos.py --check`, 2026-10-10 |
| S5 | Mesh runs: 20-36 s, peak RSS 1.52 GB, free at least 49 percent | `logs/mesh_log.tsv` |
| S6 | Distances, line nodes and quality of all 16 meshes | `logs/verify_log.tsv` |
| S7 | MOOSE input checks and the 400-way split of the largest mesh | `logs/check_decks_moose.tsv` |
| S8 | Existing two-hole mesh: band median edge 131 um at 100 um | `meshfile_2holes/mesh_2holes_3d_x1.50mm_y1.6mm_radi0.5mm.msh` |
| S9 | Frontera tree holds 0 files; no `farms-opt` anywhere on scratch | ssh listing, 2026-10-10 |
| S10 | Frontera queues close permanently on 2026-10-15 | TACC Frontera user guide, 2026-09-16 |
| S11 | Frontera MOOSE builds use gcc 9.1.0 and mvapich2-x 2.3 | `petscvariables` of `$WORK/projects/moose94b6528` |
| S12 | Sine centreline at most 0.029 mm off \eqref{eq:sine}; slices are circles | 200-slice fit of the STLs, 2026-10-10 |
| S13 | Paper setup figure: round notch tip and specimen dimensions | `Paper_ThreePointBending/Figures/problem_setup.png` |

## Appendix A: Derivations

### A.1 Throat height \eqref{eq:throat} {#a.1}

Put the contact point at the origin, with equal holes of radius $r$ centred at $(\mp r, 0)$. At height $y$ the left hole ends at $x = -r + \sqrt{r^2 - y^2}$ and the right hole starts at $x = r - \sqrt{r^2 - y^2}$. The material between them is therefore $2(r - \sqrt{r^2 - y^2})$ wide, which is zero at the contact. Requiring this width to equal $w_c$ at $y = s_0$ gives $\sqrt{r^2 - s_0^2} = r - w_c/2$, which is \eqref{eq:throat}. For $r$ = 1 mm and $w_c$ = 0.1 mm, $s_0$ = 0.312 mm. The trimmed cross-section, above and below the contact, is
$$4\int_0^{s_0} \bigl(r - \sqrt{r^2 - y^2}\bigr)\,dy = 4\Bigl[r s_0 - \tfrac{s_0}{2}\sqrt{r^2 - s_0^2} - \tfrac{r^2}{2}\arcsin(s_0/r)\Bigr] = 0.0205\ \mathrm{mm}^2,$$
or 0.16 mm$^3$ over $T$.

### A.2 Volumes \eqref{eq:volume} {#a.2}

The round notch is a rectangle $w_N \times (d_N - w_N/2)$ plus a half disk of radius $w_N/2$, which gives $A_N$ = 1.4927 mm$^2$. The flat STL notch has $A_N = w_N d_N$ = 1.6 mm$^2$, so the round tip keeps 0.107 mm$^2$ more material.

Void volumes:

- Family (c): $V_{\mathrm{void}} = 2 \pi r^2 T$.
- Family (d): $V_{\mathrm{void}} = T\bigl[\pi (r_L^2 + r_R^2) + \int_{y_0}^{y_1} \max\bigl(0,\, x_R - x_L - c_L - c_R\bigr)\,dy\bigr]$, with half chords $c_{L,R} = \sqrt{\max(0,\, r_{L,R}^2 - (y - y_c)^2)}$; the integral uses 200,000 midpoints.
- Family (e): every slice $z$ = const of \eqref{eq:sine} is a disk of area $\pi r_v^2$, so $V_{\mathrm{void}} = \pi r_v^2 T$ by Cavalieri's principle.
- Family (f) is a union instead: lower block $T(4\,\mathrm{mm}\cdot L - A_N)$, upper block $T \cdot 4\,\mathrm{mm} \cdot L$, and $n_p$ pillars of $\pi r_p^2 \cdot 0.5$ mm each.

## Appendix B: Rationale

None yet.
