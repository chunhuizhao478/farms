---
title: "Implementation Plan: SPE-HFTC acrylic-cube pulse-power cases (pure solid, 3D)"
date: 2026-10-21
geometry: margin=1in
fontsize: 11pt
colorlinks: true
header-includes: |
  \usepackage{amsmath,amssymb,booktabs}
---

## Summary

**The problem.** The existing 3D setup models a granite cylinder under ten
contiguous 5 J pulses. The SPE manuscript needs a
2-inch acrylic cube with a borehole, twenty 144 J pulses 80 s apart, a holder
clamping its ends, and the measured fracture extents after 5, 11, 15 and 20
pulses. The paper omits the energy share reaching the wall, the pressure rise
time and damping. Two existing energy postprocessors
read the wrong surfaces.

**The fix.** New mesh and inputs for the cube; paper values where given;
toughness back-calculated from the measured strength; the other missing
quantities become explicit swept parameters; a fracture-extent measure is
added; the bookkeeping is corrected.

**Expected outcome.** Each pulse costs about 700 time steps, seven tenths of
the roughly 1000 steps of the existing production run.

**Main risk / tradeoff.** Matching 69 MPa at an affordable length scale
(0.6 mm) forces a toughness of 2.7 kJ/m$^2$, eight times the literature value,
so cracks propagate too reluctantly.

**Out of scope.** Fluid coupling, field scale, the calibration campaign itself,
new C++.

## Data Flow

Eleven stages in the order the code runs. A bare name lives in memory, a path
is a file relative to `pulsepower/spe_hftc/`, and **bold** marks what ships.

**Getting ready**

1. **Build the mesh.** The gmsh script turns the cube and borehole dimensions
   of the paper's Table 2 into tetrahedra with seven named surfaces:
   `mesh/spe_cube_<tag>.msh`, plus `mesh/spe_cube_<tag>.json` with counts,
   areas and quality. [Phase\ 1](#phase-1)
2. **Derive the loading numbers.** A script solves the pulse shape rates from
   the paper's rise time and discharge duration by \eqref{eq:shape}, tabulates
   the wall pressure \eqref{eq:peak} for each energy efficiency and the
   toughness implied by the measured strength \eqref{eq:strength}: `code_3d/params_table.md`,
   `code_3d/pulse_profiles.png`. [Phase\ 2](#phase-2)

**One job (one pulse train)**

3. **Read the mesh and parameters.** `code_3d/elasticity.i` loads the `.msh`, adds
   the pinned-corner nodeset and launches `code_3d/fracture.i` on a clone of the
   same mesh. [Phase\ 3](#phase-3)
4. **Evaluate the borehole pressure.** At every quadrature point of the
   `borehole` surface the pulse function gives $p(t,\mathbf{x})$ by
   \eqref{eq:pulse} and \eqref{eq:peak}: a traction on the wall and the aux
   field `pulse_load_aux`. [Phase\ 3](#phase-3)
5. **Advance the solid.** One Newmark step of \eqref{eq:momentum} with the
   degraded stress \eqref{eq:stress} and Rayleigh damping, side faces free and
   ends held vertically: `disp_*`, `vel_*`, `accel_*`, `psie_active`.
   [Phase\ 3](#phase-3)
6. **Advance the damage.** The sub-app solves \eqref{eq:pff} under the
   irreversibility bounds, driven by `psie_active`: `d` and
   `dissipated_energy_density`, copied back. [Phase\ 3](#phase-3)
7. **Iterate and step.** Stages 5 and 6 repeat until the fixed-point
   tolerance; the step size follows the Newton iteration count, is capped
   during each rise, lands exactly on every pulse onset and restarts at
   $\Delta t_0$ there. [Phase\ 3](#phase-3)
8. **Measure.** Fracture extent \eqref{eq:extent}, crack area \eqref{eq:area}
   and the energy budget \eqref{eq:budget} every step:
   **`code_3d/cases/<case>/elasticity_csv.csv`**; damage snapshots
   **`code_3d/cases/<case>/elasticity_out.e`**. [Phase\ 3](#phase-3)

**Across jobs and after**

9. **Check a small run.** The smoke case on the coarse mesh runs one pulse,
   undamped and damped, and compares the wall pressure, the energy budget,
   the sideset areas and the extent metric with their expected values:
   report `code_3d/smoke/SMOKE_REPORT.md`. [Phase\ 4](#phase-4)
10. **Submit, continue and fetch.** SLURM scripts split the mesh, run, continue
    with `--recover` past each allocation limit, and rsync results as the
    existing cases do. [Phase\ 5](#phase-5)
11. **Compare with the experiment.** A script reads each case's CSV, takes the
    extent at the end of each completed pulse window and plots it against the
    four measured points: **`analysis/fracture_extent.png`**,
    **`analysis/energy_budget.png`**. [Phase\ 5](#phase-5)

Stages 3 to 8 are the existing pipeline with new parameters and two added
measures; stages 1, 2, 9 and 11 are new.

## Notation

| Symbol | Meaning | Units | Code variable |
| ------------- | --------------------------------------------- | ----------- | --------------------------- |
| $\rho$ | mass density | kg/m$^3$ | `density` |
| $E$, $\nu$ | Young's modulus, Poisson ratio | Pa, -- | `E`, `nu` |
| $K$, $G$, $\lambda$ | bulk, shear and first Lamé moduli | Pa | `K`, `G`, `lambda` (in material) |
| $c_p$, $c_s$ | P- and S-wave speeds | m/s | `Cp`, `Cs` |
| $f_t$ | measured tensile strength (paper) | Pa | `ft` (comment only) |
| $G_c$ | fracture toughness, back-calculated from $f_t$ and $l$ | J/m$^2$ | `Gc_const` |
| $l$ | phase-field regularisation length | m | `l` |
| $l_{ch}$ | Irwin length $EG_c/f_t^2$ | m | `l_ch` (script) |
| $d$ | damage field, 0 intact, 1 broken | -- | `d` |
| $\alpha(d)$, $c_0$ | crack geometric function $\alpha = d$ and its normalisation $8/3$ | -- | `alpha`, `c0` |
| $g(d)$, $\eta_r$ | degradation $(1-d)^2(1-\eta_r)+\eta_r$ and residual stiffness | -- | `g`, `eta` |
| $\sigma_c$ | AT1 homogeneous tensile strength | Pa | `sigma_c` (script) |
| $\mathbf{u}$, $\mathbf{v}$, $\mathbf{a}$, $u_z$ | displacement, velocity, acceleration; vertical displacement | m, m/s, m/s$^2$, m | `disp_*`, `vel_*`, `accel_*`, `disp_z` |
| $\boldsymbol\varepsilon$, $\boldsymbol\sigma$ | small strain, Cauchy stress | --, Pa | `mechanical_strain`, `stress` |
| $\boldsymbol\sigma^{\pm}$, $\boldsymbol\varepsilon^{+}$ | spectral tensile/compressive stress parts, tensile strain part | Pa, -- | internal to `elasticity` |
| $\psi^+$, $\psi_e$ | active and total elastic energy density | J/m$^3$ | `psie_active`, `psie` |
| $\psi_d$ | crack (dissipated) energy density | J/m$^3$ | `dissipated_energy_density` |
| $\beta_N$, $\gamma_N$ | Newmark parameters | -- | `newmark_beta`, `newmark_gamma` |
| $\zeta$ | stiffness-proportional Rayleigh damping coefficient | s | `zeta` |
| $\xi$, $\omega$ | damping ratio, angular frequency | --, rad/s | script only |
| $p(t,\mathbf{x})$ | borehole-wall pressure | Pa | `func_pulse` |
| $P_p$ | peak pressure of the empirical law | Pa | `Pp` (script) |
| $a$, $b$ | pulse shape rates | 1/s | `shape_alpha`, `shape_beta` |
| $t_r$, $t_d$ | pressure rise time, discharge duration to 10 % | s | `t_r`, `t_dur` |
| $\tau$ | time inside the current pulse window | s | `mod_time` |
| $T_w$, $N_p$ | pulse window, pulses per job | s, -- | `Tw`, `Np` |
| $E_M$ | stored energy per pulse | kJ | `EM` |
| $\eta$ | fraction of $E_M$ converted to borehole pressure | -- | `eta` |
| $B$, $\kappa$, $m$ | constants of the peak-pressure law | -- | `base_factor`, `kappa`, `m_exp` |
| $D(\mathbf{x})$ | distance from $\mathbf{x}$ to the discharge segment | mm | `Dsensor` |
| $\mathbf{x}_c$, $L_g$ | discharge centre, electrode gap length | m | `discharge_center`, `gap` |
| $L_x$, $L_y$, $L_z$ | specimen dimensions | m | `Lx`, `Ly`, `Lz` |
| $r_b$ | borehole radius | m | `rb` |
| $h_f$, $h_c$, $r_f$ | fine and coarse element size, fine-zone radius | m | `--h-fine`, `--h-coarse`, `--r-fine` |
| $N_e$ | number of tetrahedra | -- | `n_tets` (json) |
| $d_{th}$ | damage threshold that counts as fractured | -- | `d_th` |
| $L_f$ | fracture extent, longest crack from the wall in top view | m | `fracture_extent` |
| $V_d$ | damaged volume | m$^3$ | `damaged_volume` |
| $A_f$ | crack-area estimate | m$^2$ | `crack_area` |
| $W_{ext}$ | work of the pressure forces | J | `external_work` |
| $E_{in}$, $E_{tot}$ | energy input, energy stored plus dissipated | J | `full_input_energy`, `full_energy` |
| $\mathcal{K}$, $\Psi_e$, $\Psi_d$ | kinetic, elastic and dissipated energy | J | `solid_kinetic_energy_total`, `solid_elastic_energy_total`, `dissipated_energy_total` |
| $\Delta t_0$, $\Delta t_r$, $\Delta t_{max}$ | onset, rise-phase cap and maximum time step | s | `dt0`, `dt_rise`, `dt_max` |
| $\mathbf{n}$ | outward unit normal | -- | `_normals` |
| $H$ | Heaviside step | -- | `if(...)` in `ParsedAux` |
| $\epsilon$ | tolerance of the energy-budget check | -- | `check_smoke.py` |
| $x$, $y$, $z$, $r$ | coordinates (borehole axis along $z$) and radius $\sqrt{x^2+y^2}$ | m | `x`, `y`, `z` in `ParsedAux` |
| $k$ | pulse index | -- | `pulse` (script) |
| $w$, $\mathbf{f}$ | test function, nodal force vector | --, N | `_test`, `fx fy fz` |
| $\Delta t$ | current time step | s | `_dt` |
| $d_{\text{old}}$ | damage at the previous step (lower bound) | -- | `initial_damage_aux`, VI lower bound |
| $\langle\cdot\rangle_+$ | Macaulay bracket, $\max(\cdot,0)$ | -- | `RaccoonUtils::Macaulay` |
| $\Pi$ | regularised total energy functional | J | none (Appendix A.1) |
| $\varepsilon$, $s$ | scalar uniaxial strain; integration variable | -- | none (Appendix A) |
| $f(\tau)$, $\tau_p$ | unnormalised pulse shape and its peak time | --, s | `tpeak` (script) |
| $k_m$, $m_0$, $c$, $\mathbf{K}$ | modal stiffness, mass, damping constant; stiffness matrix | N/m, kg, N s/m, N/m | none (Appendix A.4) |
| $\sigma_{\theta\theta}$ | hoop stress | Pa | none (Appendix A.8) |
| $\mathbb{C}$, $\boldsymbol\sigma_v$ | degraded elasticity tensor, damping stress $\zeta\dot{\boldsymbol\sigma}$ | Pa, Pa | internal to the kernel |
| $\mathbf{M}$, $\mathbf{C}_d$, $q$ | mass matrix, damping matrix, modal coordinate | kg, N s/m, m | none (Appendix B) |
| $\tan\delta$ | loss factor of a Kelvin--Voigt solid, $\zeta\omega$ | -- | none (Appendix B) |

## Formulation

\begin{equation}\label{eq:momentum}
\rho\,\ddot{\mathbf{u}} = \nabla\cdot\left(\boldsymbol\sigma + \zeta\,\dot{\boldsymbol\sigma}\right)
\ \text{in } \Omega,\qquad
\boldsymbol\sigma\mathbf{n} = -p\,\mathbf{n}\ \text{on the borehole},\qquad
u_z = 0,\ (\boldsymbol\sigma\mathbf{n})_{x,y} = 0\ \text{on the ends},\qquad
\boldsymbol\sigma\mathbf{n} = \mathbf{0}\ \text{on the side faces}
\end{equation}
*From [M1]; $\zeta>0$ is a modeling assumption standing in for the material
damping of acrylic: with $\zeta = 0$ the cube rings forever. The end condition
is the holder of Fig. 2a as decided by the user.*

\begin{equation}\label{eq:stress}
\boldsymbol\sigma = g(d)\,\boldsymbol\sigma^{+}(\boldsymbol\varepsilon) + \boldsymbol\sigma^{-}(\boldsymbol\varepsilon),\qquad
\psi^+ = \tfrac12\lambda\langle\operatorname{tr}\boldsymbol\varepsilon\rangle_+^2 + G\,\boldsymbol\varepsilon^{+}\!:\!\boldsymbol\varepsilon^{+}
\end{equation}
*From [S1], spectral branch.*

\begin{equation}\label{eq:pff}
\frac{G_c}{c_0}\left(\frac{1}{l} - 2l\,\Delta d\right) + g'(d)\,\psi^+ = 0,\qquad
c_0 = \tfrac{8}{3},\qquad d_{\text{old}} \le d \le 1
\end{equation}
*Derived* -- Appendix A.1.

\begin{equation}\label{eq:strength}
\sigma_c = \sqrt{\frac{3\,E\,G_c}{8\,l}}
\quad\Longleftrightarrow\quad
G_c = \frac{8\,l\,\sigma_c^2}{3\,E},\qquad \sigma_c := f_t
\end{equation}
*Derived* -- Appendix A.2. The plan fixes $\sigma_c = f_t$ and takes $G_c$ from
the right-hand form.

\begin{equation}\label{eq:pulse}
p(t,\mathbf{x}) = P_p(\mathbf{x})\,\frac{e^{-a\tau}-e^{-b\tau}}{e^{-a t_r}-e^{-b t_r}},\qquad
\tau = t \bmod T_w,\qquad t \le N_p T_w
\end{equation}
*From [S2].*

\begin{equation}\label{eq:peak}
P_p(\mathbf{x}) = 10^{5}\,B\,\frac{(\eta\,E_M)^{\kappa}}{D(\mathbf{x})^{m}}
\end{equation}
*From [S2]; $B=8000$, $\kappa=0.35$, $m=0.25$ are the user's earlier calibration (Open Questions).*

\begin{equation}\label{eq:shape}
\frac{\ln(b/a)}{b-a} = t_r,\qquad \frac{p(t_d)}{P_p} = 0.1
\end{equation}
*Definition* -- peak at the rise time, 10 % left at the discharge duration; if
the pressure history differs from the current history the impulse per pulse is
wrong by the same factor. Appendix A.3.

\begin{equation}\label{eq:extent}
L_f = \max_{\text{nodes with } d \ge d_{th}} \sqrt{x^2+y^2}\;-\;r_b,\qquad
V_d = \int_\Omega H(d-d_{th})\,dV
\end{equation}
*Definition* -- the paper's "longest fracture from borehole boundary to tip",
in top view.

\begin{equation}\label{eq:area}
A_f = \frac{\Psi_d}{G_c},\qquad
\Psi_d = \int_\Omega \frac{G_c}{c_0}\left(\frac{\alpha(d)}{l} + l\,|\nabla d|^2\right) dV
\end{equation}
*From [S8] for $\Psi_d$; $A_f$ is a definition* -- exact for a developed AT1 profile.

\begin{equation}\label{eq:budget}
E_{in} = -W_{ext},\qquad E_{tot} = \mathcal{K}+\Psi_e+\Psi_d,\qquad
|E_{tot}-E_{in}| \le \epsilon\,\max_t E_{in}\ \ (\zeta = 0)
\end{equation}
*From [S4] for the work term; the clamped ends do no work; the inequality is
the acceptance check.* With $\zeta>0$ the gap is the damped energy.

\begin{equation}\label{eq:rayleigh}
\xi(\omega) = \tfrac12\,\zeta\,\omega
\end{equation}
*Derived* -- Appendix A.4.

- First-order Lagrange tetrahedra.
- Newmark average acceleration, $\beta_N = 1/4$, $\gamma_N = 1/2$.
- Damage as a variational inequality; fields coupled by fixed-point iteration.

## Approach

Copy the `3d_solid/em_5J` inputs to `spe_hftc/code`, then replace the mesh, the
material constants, the boundary names and the pulse parameters, and add
Rayleigh damping, a pulse-onset sync function plus a rise-phase step cap for
the timestepper, the extent and area postprocessors, and the corrected
boundaries of the two energy postprocessors. The holder plates become
$u_z = 0$ on the end faces. A gmsh
Python script `mesh/make_spe_cube_mesh.py` builds the cube with named surfaces
and a cylindrical refinement zone; `mesh/make_spe_section_2d.py` builds the
plane-strain section at the gap height, and `code_2d/` holds the matching 2D
inputs. The gotchas of the existing setup are recorded in [E1].

## Called Code

\begingroup\footnotesize

| Function | What it computes | Governing equation | Verified how |
| ------------------------------------------ | -------------------------------- | -------- | ----------------- |
| `ElkPulseLoadExperiment::value` | $p(t,\mathbf{x})$ with $D$ to the gap segment in mm | \eqref{eq:pulse}, \eqref{eq:peak} | read source [S2] |
| `Pressure` action | traction $-p\mathbf{n}$ on `borehole`; nodal forces into `fx fy fz` | \eqref{eq:momentum} | read `validParams`; production use [S7] |
| `ADDynamicStressDivergenceTensors` | $\nabla w\cdot[(1+\zeta/\Delta t)\boldsymbol\sigma-(\zeta/\Delta t)\boldsymbol\sigma_{\text{old}}]$ at zero HHT | \eqref{eq:momentum} | read source [M1] |
| `ADInertialForce` | $w\,\rho\,\mathbf{a}$ from the Newmark update | \eqref{eq:momentum} | read source [M1] |
| `NewmarkBeta`, `NewmarkAccelAux`, `NewmarkVelAux` | Newmark update of $\mathbf{u},\mathbf{v},\mathbf{a}$ | \eqref{eq:momentum} | production use [S7] |
| `SmallDeformationIsotropicElasticity` (SPECTRAL) | $\boldsymbol\sigma$, $\psi^+$, $\psi_e$ | \eqref{eq:stress} | read source [S1] |
| `ComputeSmallDeformationStress`, `ADComputeSmallStrain` | calls the model; $\boldsymbol\varepsilon = \operatorname{sym}\nabla\mathbf{u}$ | \eqref{eq:stress} | read source [S1]; production use |
| `PowerDegradationFunction` | $g(d)$, $g'(d)$ | \eqref{eq:stress}, \eqref{eq:pff} | read source [S8] |
| `CrackGeometricFunction` | $\alpha(d)=d$ and $c_0$ by adaptive trapezoid | \eqref{eq:pff} | read source [S8] |
| `ADPFFDiffusion`, `ADPFFSource` | $(\nabla w,\,2G_c l/c_0\,\nabla d)$ and $(w,\,\partial\psi/\partial d)$ | \eqref{eq:pff} | read source [S8] |
| `ADDerivativeParsedMaterial` (`psi`) | $\psi = \alpha G_c/(c_0 l) + g\,\psi^+$ and $\partial\psi/\partial d$ | \eqref{eq:pff} | production use [S7] |
| `VariableConstantIrreversibleBounds`, `ConstantBounds` | $d_{\text{old}} \le d \le 1$ for the VI solver | \eqref{eq:pff} | read source [S8] |
| `CrackDissipatedEnergyDensity` | $\psi_d$ | \eqref{eq:area} | read source [S8] |
| `FirstStepElementIntegralVariablePostprocessor` | freezes $\int\psi_d$ of the first step | \eqref{eq:budget} | read source [S4] |
| `FarmsExternalWork` | $\sum_{\text{boundary nodes}}\mathbf{f}\cdot\dot{\mathbf{u}}\,\Delta t$, accumulated | \eqref{eq:budget} | read source [S4] |
| `ADKineticEnergyAux`, `ElementIntegralVariablePostprocessor` | $\tfrac12\rho|\mathbf{v}|^2$ and its integral | \eqref{eq:budget} | read source [M1] |
| `ADElementIntegralMaterialProperty` | $\Psi_e = \int\psi_e\,dV$ | \eqref{eq:budget} | production use [S7] |
| `DirichletBC` | $u_z = 0$ on the ends; $u_x = u_y = 0$ at the corner node | \eqref{eq:momentum} | production use [S7] |
| `GenericFunctionMaterial`, `PiecewiseConstant` (`zeta_fn`) | the `zeta` property as a function of time (quench schedule) | \eqref{eq:momentum} | 2D run of 2026-10-21 (`--check-input` and 10 $\mu$s) |
| `FarmsIterationAdaptiveDT` | iteration-adaptive $\Delta t$, cap, onset sync, post-sync $\Delta t_0$, postprocessor cap | none -- stepping | read source [S5]; this combination **NOT verified** by a run |
| `PiecewiseLinear`, `PiecewiseConstant`, `FunctionValuePostprocessor` | onset list; rise-phase $\Delta t_r$ schedule | none -- stepping | read source (`direction` semantics); **NOT verified** by a run |
| `ParsedAux` (`use_xyzt`), `NodalExtremeValue`, `ElementExtremeValue`, `ParsedPostprocessor`, `AreaPostprocessor` | $L_f$, $V_d$, $A_f$, wall-pressure peak, sideset area | \eqref{eq:extent}, \eqref{eq:area} | classes exist in the framework; **NOT verified** by a run |
| `FileMeshGenerator`, `ExtraNodesetGenerator` | mesh and corner nodeset | none | production use [S7] |
| `MultiAppCopyTransfer`, `MultiAppPostprocessorTransfer`, `Receiver` | field and scalar exchange between apps | none | production use [S7] |
| `MeshSize`, `ErrorPsiMeasure` | diagnostics kept for parity | none | read source [S8] |
| gmsh API: `occ.addBox`, `occ.addCylinder`, `occ.cut`, `mesh.field`, `getCenterOfMass`, `getMass` | geometry, size field, face classification, areas | none | gmsh 4.15.1 importable; script **NOT verified** |

\endgroup

## Constraints

- No new C++: only objects registered in the `farms-opt` build on Frontera (`farms_cdms_04192026`) and locally.
- Main app and sub-app load the identical mesh (`clone_parent_mesh`, `MultiAppCopyTransfer`).
- Production mesh at most 20 M tetrahedra; the existing production ran 16.4 M on 500 cores for 48 h.
- $\Delta t \le \Delta t_r = 3\times10^{-8}$ s during the first 0.6 $\mu$s of each pulse; $\Delta t_{max} = 10^{-7}$ s.
- $l \ge 2h_f$ inside the fine zone.
- Local runs need the `mpicxx` shim in `tools/jit_shim` on `PATH` so the FParser JIT links (macOS toolchain, [T1]); AD parsed materials cannot run without the JIT.

## Phases at a Glance

\begingroup\small

| # | After this phase | Stages | Files touched | Done when |
| --- | ------------------------------ | ------ | ------------------------------ | ---------------------------- |
| [1](#phase-1) | A script makes cube meshes with named faces at any resolution | 1 | `mesh/make_spe_cube_mesh.py`, `mesh/README.md`, `mesh/*.msh`, `mesh/*.json` | smoke and production meshes exist; counts within 20 % of Appendix A.5 |
| [2](#phase-2) | The loading and strength numbers are tabulated and plotted | 2 | `code_3d/derive_params.py`, `code_3d/params_table.md`, `code_3d/pulse_profiles.png` | table matches Appendix A to 0.1 % |
| [3](#phase-3) | The acrylic-cube inputs parse and load the mesh | 3--8 | `code_3d/elasticity.i`, `code_3d/fracture.i` | `--check-input` and `--mesh-only` pass |
| [4](#phase-4) | A one-pulse coarse run balances energy and reports the extent | 9 | `code_3d/smoke/run_smoke.sh`, `code_3d/smoke/check_smoke.py`, `code_3d/smoke/SMOKE_REPORT.md` | budget error $\le 5$ %; pressure peak within 2 % |
| [5](#phase-5) | Cases run on Frontera as continuous pulse trains and plot against the paper | 10--11 | `code_3d/make_cases.py`, `code_3d/case_template/*.sbatch`, `code_3d/submit_all.sh`, `code_3d/sync_results.sh`, `code_3d/post_fracture_extent.py` | recover test passes; plot shows the four paper points |

\endgroup

## Open Questions

| Question | What it blocks | How to resolve |
| --------------------------------------------- | ------------------ | ---------------------------------- |
| How to bring the cube to rest between pulses: constant Rayleigh $\zeta = 2\times10^{-7}$ s also damps the 0.3 $\mu$s rise (damping stress $\zeta/t_r = 0.67$ of the elastic stress), whereas a quench ($\zeta = 10^{-5}$ s switched on 25 $\mu$s after each onset, $T_w = 60\ \mu$s) leaves the pulse untouched; which one? | [Phase\ 3](#phase-3) damping and window | user decision; recommendation is the quench, see Appendix B |
| Do the new metric postprocessors and the step schedule behave as specified? | [Phase\ 4](#phase-4) | the smoke run checks (a)--(e) |
| Does the production CSV of `3d_solid/em_5J` show `full_input_energy` near zero (confirming the boundary swap, [S4])? | interpretation of existing cmame results only | inspect the CSV on Frontera or the external drive |
| What is the wall time per step on Frontera for the 16.4 M mesh? | sizing of [Phase\ 5](#phase-5) | read the `.o` logs of the `em_5J` job on Frontera |

## Phase 1: Mesh generator {#phase-1}

**In one sentence.** Turns the specimen dimensions into a tetrahedral mesh of
the cube minus the borehole, fine over the whole height inside a cylinder
around the borehole and coarse in the corners, with seven named surfaces.

**Implements [Data Flow](#data-flow) stage(s):** 1

**Files.** `mesh/make_spe_cube_mesh.py` (new) -- the generator.
`mesh/README.md` (new) -- the two command lines. `mesh/spe_cube_smoke.msh`,
`mesh/spe_cube_h0p3.msh` and their `.json` (generated).

**Requirements.**

1. `def build(lx: float, ly: float, lz: float, rb: float, h_fine: float, h_coarse: float, r_fine: float, out: str, msh_version: str = "4.1", threads: int = 8) -> dict`
   writes `<out>.msh` (format 4.1, as the production mesh) and `<out>.json`,
   and returns the json dict. CLI flags mirror the arguments
   (`--lx 51.11e-3 --ly 51.92e-3 --lz 49.14e-3 --rb 4.15e-3 --h-fine 0.3e-3 --h-coarse 1.0e-3 --r-fine 26e-3 --out spe_cube_h0p3`).
2. Geometry with `gmsh.model.occ`: box $[-L_x/2, L_x/2]\times[-L_y/2, L_y/2]\times[0, L_z]$
   minus the cylinder of radius $r_b$ along $z$ through $(0,0)$, by `occ.cut`
   then `occ.synchronize`.
3. Physical groups by geometry, never by OCC tag: classify each surface by
   `getCenterOfMass` (face centre within $10^{-6}$ m of $\pm L_x/2$, $\pm L_y/2$,
   $0$ or $L_z$) into `xmin xmax ymin ymax zmin zmax`; the remaining surface,
   whose bounding box has $x$-extent $2r_b$, is `borehole`; the volume is
   `acrylic`. Raise if the count is not 7 + 1.
4. Size field: a `MathEval` background field
   $h_f + (h_c - h_f)\,\max(0,\min(1,(\sqrt{x^2+y^2}-r_f)/3\,\text{mm}))$;
   set `Mesh.MeshSizeExtendFromBoundary = 0`, `Mesh.MeshSizeFromPoints = 0`,
   `Mesh.MeshSizeFromCurvature = 0`; `Mesh.Algorithm3D = 10` (HXT) with
   `General.NumThreads = threads`; `Mesh.Optimize = 1`.
5. The json records the CLI arguments, `n_nodes`, `n_tets`, triangles per
   surface group, each group's area from `occ.getMass`, the volume, and the
   min and mean `gamma` quality from `gmsh.model.mesh.getElementQualities`.
6. Two configurations documented in `mesh/README.md`: `spe_cube_smoke`
   ($h_f = h_c = 1.5$ mm) and `spe_cube_h0p3` ($h_f = 0.3$ mm, $h_c = 1.0$ mm,
   $r_f = 26$ mm). Optional `spe_cube_h0p4` ($h_f = 0.4$ mm) for sweeps.

**Interfaces.** Both inputs of [Phase\ 3](#phase-3) read the `.msh` through
`FileMeshGenerator` and use the surface and block names.

**Edge cases.**

- `occ.cut` returns a face count other than 7 $\to$ raise, listing the face centres.
- $h_f > r_b/8$ $\to$ warn: fewer than 25 segments around the borehole.
- Minimum `gamma` below 0.1 $\to$ rerun with `Mesh.OptimizeNetgen = 1` on the smoke mesh only (too slow for 17 M).

**Acceptance criteria.**

- [ ] Both meshes and json files exist; $N_e$ within 20 % of Appendix A.5 (0.17 M and 17 M).
- [ ] Exactly 7 surface groups; `borehole` area within 0.5 % of $2\pi r_b L_z = 1.2813\times10^{-3}$ m$^2$; volume within 0.5 % of $1.2774\times10^{-4}$ m$^3$.
- [ ] Minimum `gamma` $\ge 0.1$ on both meshes.
- [ ] `gmsh.open` reloads each file without error.

**Depends on:** nothing

## Phase 2: Loading and property derivation {#phase-2}

**In one sentence.** Converts the paper's pulse and material data into the
numbers the inputs need and writes them as a table and a plot.

**Implements [Data Flow](#data-flow) stage(s):** 2

**Files.** `code_3d/derive_params.py` (new) -- the derivation; `code_3d/params_table.md`
and `code_3d/pulse_profiles.png` (generated).

The quantities the model needs, and where each stands against the paper [P1]:

| Quantity | Symbol | Paper value | Status and plan default |
| -------------------- | ------------ | ------------------------- | ------------------------------------------- |
| density | $\rho$ | 1190 kg/m$^3$ (Table 1) | given |
| Young's modulus | $E$ | 2.80 GPa | given |
| Poisson ratio | $\nu$ | 0.35 | given |
| tensile strength | $f_t$ | 69 MPa | given; sets $G_c$ through \eqref{eq:strength} |
| compressive strength | -- | 124 MPa | given; unused (no compressive failure in the model) |
| fracture toughness | $G_c$ | none | **back-calculated** by \eqref{eq:strength}: 2721 J/m$^2$ at $l = 0.6$ mm; literature 350 J/m$^2$ [HM12] kept as a variant |
| regularisation length | $l$ | none (numerical) | 0.6 mm, set by the mesh; $\sigma_c = f_t$ by construction |
| specimen size | $L_x, L_y, L_z$ | 51.11, 51.92, 49.14 mm (Table 2) | given |
| borehole radius | $r_b$ | 4.15 mm (8.3 mm, Fig. 2 caption; Table 2 says 7.94 mm) | 8.3 mm per the user (2026-10-21) |
| electrode gap | $L_g$ | 7.17 mm | given |
| gap position | $\mathbf{x}_c$ | none | $(0,0,L_z/3)$ from the bottom face, per the user (2026-10-21) |
| stored energy | $E_M$ | 144 J (60 kV, 80 nF) | given |
| conversion efficiency | $\eta$ | none | free parameter, swept 0.1--1 to see how the pattern changes, per the user |
| peak-pressure constants | $B,\kappa,m$ | none | from the existing setup (8000, 0.35, 0.25), kept for now per the user (2026-10-21) |
| pressure rise time | $t_r$ | current rise 230--350 ns | 0.3 $\mu$s, equal to the current rise, per the user |
| discharge duration | $t_d$ | 3--12 $\mu$s to 10 % | given as a range; 6 $\mu$s, sweep 3 and 12 |
| pulse interval | -- | 80 s | replaced by $T_w = 60\ \mu$s with a damping quench after 25 $\mu$s |
| pulses | $N_p$ | 20, in 4 batches of 5 | one continuous train, no reset between pulses; the count is not binding, per the user |
| material damping | $\zeta$ | none | **missing**; $2\times10^{-7}$ s, \eqref{eq:rayleigh} |
| holder constraint | -- | plates on top and bottom (Fig. 2a) | ends held vertically, $u_z = 0$ on `zmin` and `zmax`; side faces free, per the user |
| residual stiffness, Newmark | $\eta_r$, $\beta_N$, $\gamma_N$ | none | $10^{-6}$, 1/4, 1/2 as in the existing setup |
| borehole fluid, electrodes | -- | 2 % NaCl, Cu-101 | not needed for a pure solid |

**Requirements.**

1. `def shape_params(t_r: float, t_d: float, level: float = 0.1) -> tuple[float, float]`
   solves \eqref{eq:shape}: inner `brentq` for $b \in (a(1+10^{-9}), 10^{10})$
   from the peak-time identity at fixed $a$; outer `brentq` for
   $a \in (10^2, 0.999/t_r)$ from the decay condition. Raise `ValueError` when
   no bracket exists.
2. `def peak_pressure(EM_kJ: float, D_mm: float, eta: float, B: float = 8000.0, kappa: float = 0.35, m: float = 0.25) -> float`
   implements \eqref{eq:peak}; `def distance_to_gap(x, y, z, xc, gap) -> float`
   reproduces $D$ of [S2] (radial inside the gap's $z$-range, distance to the
   nearer gap end outside it), in mm.
3. `def at1_strength(E, Gc, l) -> float` and its inverse
   `def gc_from_strength(E, ft, l) -> float`, both \eqref{eq:strength};
   `def wave_speeds(E, nu, rho) -> tuple[float, float]`;
   `def irwin_length(E, Gc, ft) -> float`; `def rayleigh_zeta(xi, f_hz) -> float` from \eqref{eq:rayleigh}.
4. `def pulse_profile(t, a, b, t_r) -> ndarray` and a figure of the three
   durations for $t_r = 0.3\ \mu$s over 0--40 $\mu$s.
5. `params_table.md` contains: $K$, $G$, $c_p$, $c_s$; $G_c$ from $f_t = 69$ MPa
   for $l \in \{0.3, 0.5, 0.6, 0.8\}$ mm with its ratio to 350 J/m$^2$;
   $\sigma_c$ at 350 J/m$^2$ for the same $l$ (the literature-toughness
   variant); $(a, b)$ for $t_d \in \{3, 6, 12\}\ \mu$s; $P_p$ at
   the gap ($D = r_b$) and at the near and far borehole ends ($D = 13.45$ and 29.47 mm) for
   $\eta \in \{0.1, 0.3, 0.5, 1\}$; the step budget per pulse of Appendix A.7.
6. Self-check on the existing parameters ($a = 4.658\times10^5$,
   $b = 4.661\times10^5$, $t_r = 3\ \mu$s; $E_M = 0.005$ kJ, $D = 1.6$ mm):
   peak at 2.15 $\mu$s with $1.065\,P_p$, 10 % at 10.5 $\mu$s, $P_p = 111$ MPa.

**Interfaces.** The numbers are copied by hand into the header of
`code_3d/elasticity.i` in [Phase\ 3](#phase-3).

**Edge cases.**

- $t_r \ge t_d$ $\to$ `ValueError`.
- $\eta E_M \le 0$ or $D = 0$ $\to$ `ValueError`.

**Acceptance criteria.**

- [ ] $(a, b)$ reproduce Table A.3, $G_c(l)$ reproduces Table A.2 and $P_p$ reproduces Table A.6, each within 0.1 %.
- [ ] The self-check prints the four production reference numbers of requirement 6.

**Depends on:** nothing

## Phase 3: Input files for the acrylic cube {#phase-3}

**In one sentence.** Runs a train of $N_p$ pulses on the cube: reads the mesh,
applies the pulse pressure to the borehole wall, advances the damped solid and
the damage field, and writes energies, extent and damage snapshots.

**Implements [Data Flow](#data-flow) stage(s):** 3--8

**Files.** `code_3d/elasticity.i` (new, derived from `3d_solid/em_5J/elasticity.i`)
-- main app. `code_3d/fracture.i` (new, from `3d_solid/em_5J/fracture.i`) -- sub-app.
`code_2d/elasticity.i`, `code_2d/fracture.i` (new) -- the plane-strain section
twin: two displacement components, pressure on the `borehole` curve, discharge
centre at the origin, no end constraint; same parameters otherwise.
`code_3d/fracture_weibull.i`, `code_2d/fracture_weibull.i` (new) -- sub-app
variant with an element-wise Weibull-scattered $G_c$ (shape 10, scale set so
the mean is `Gc_const`, fixed seed) that breaks the axial symmetry so radial
cracks can localize; `fracture_gaussian.i` (new, both folders) -- the same with
a truncated normal (12 % standard deviation, bounds 0.3 and 2 times the mean).
Both are selected with `MultiApps/fracture/input_files=<variant>`. The base
`fracture.i` carries the same toughness field with a constant value, and the
main app copies it to the output together with the particle velocity
components, their magnitude and the maximum particle velocity `v_max`.
`code_2d/cases/discovery_gaussian_eta1/` (new) -- SLURM script and README for
the first cluster run on Discovery: Gaussian scatter, $\eta = 1$, 20 pulses.

**Requirements.**

1. Header variables, all overridable from the command line:
   `E = 2.80e9`, `nu = 0.35`, `density = 1190`, `ft = 69e6`, `l = 6e-4`,
   `Gc_const = ${fparse 8*l*ft*ft/(3*E)}` (2721 J/m$^2$; comment: eight times
   the literature toughness, Appendix B);
   `Lx = 51.11e-3`, `Ly = 51.92e-3`, `Lz = 49.14e-3`, `rb = 4.15e-3`,
   `zc = ${fparse Lz/3}`, `gap = 7.17e-3`; `EM = 0.144`, `eta = 0.3`,
   `base_factor = 8000`, `kappa = 0.35`, `m_exp = 0.25`, `t_r = 3e-7`,
   `shape_alpha = 4.1032e5`, `shape_beta = 1.1529e7`, `Tw = 6e-5`, `Np = 20`;
   `t_free = 2.5e-5`, `zeta_pulse = 0`, `zeta_quench = 1e-5` (the quench;
   constant damping is `zeta_pulse = zeta_quench`, none is both zero),
   `hht_alpha = 0`; `dt0 = 1e-8`, `dt_rise = 3e-8`, `t_rise_hold = 6e-7`,
   `dt_max = 1e-7`; `d_th = 0.5`; `mesh_file = ../mesh/spe_cube_h0p3.msh`.
2. Mesh: `FileMeshGenerator` with `file = ${mesh_file}`; `ExtraNodesetGenerator`
   with `coord = '${fparse -Lx/2} ${fparse -Ly/2} 0'`, `use_closest_node = true`,
   `new_boundary = corner_ptr`; no subdomain box generators. The block in
   `fracture.i` is identical.
3. Materials as in the existing input with the new constants; spectral
   elasticity; power degradation with $\eta_r = 10^{-6}$; `zeta_mat` of type
   `GenericFunctionMaterial` with `prop_names = zeta`, `prop_values = zeta_fn`.
4. Kernels: the three `ADDynamicStressDivergenceTensors` get `zeta = zeta`
   (the property of requirement 3) and `alpha = ${hht_alpha}`;
   `ADInertialForce` unchanged.
5. Functions: `func_pulse` of type `ElkPulseLoadExperiment` with
   `shape_param_alpha = ${shape_alpha}`, `shape_param_beta = ${shape_beta}`,
   `rise_time = ${t_r}`, `single_pulse_duration = ${Tw}`, `EM = ${EM}`,
   `gap = ${gap}`, `convert_efficiency = ${eta}`, `fitting_param_alpha = ${kappa}`,
   `fitting_param_exponent = ${m_exp}`, `base_factor = ${base_factor}`,
   `discharge_center = '0 0 ${zc}'`, `number_of_pulses = ${Np}`;
   `pulse_sync` of type `PiecewiseLinear` with `x` = the onsets $kT_w$,
   $k = 0..40$ and `y` = $0..40$; `dt_limit` of type `PiecewiseConstant`
   (`direction = LEFT`) with `x` = $\{kT_w,\ kT_w + 0.6\ \mu\text{s}\}_{k=0..40}$ and
   `y` alternating `${dt_rise}`, `${dt_max}`; `zeta_fn` of type
   `PiecewiseConstant` (`direction = LEFT`) with `x` = the onsets and the
   onsets plus `t_free`, `y` alternating `${zeta_pulse}`, `${zeta_quench}`;
   points past `end_time` are never reached.
6. BCs: `[Pressure/pressure_borehole]` with `boundary = borehole`,
   `function = func_pulse`, `save_in_disp_x/y/z = fx/fy/fz`; `fix_z_ends` of
   type `DirichletBC` on `disp_z` with `boundary = 'zmin zmax'`, `value = 0`
   (the holder plates); `fix_corner_x` and `fix_corner_y` Dirichlet 0 on
   `corner_ptr` (in-plane rigid modes). No dashpots; the side faces are free.
7. Executioner as existing, with `end_time = ${fparse Np*Tw}`;
   `[TimeStepper]` of type `FarmsIterationAdaptiveDT` with `dt = ${dt0}`,
   `max_time_step_bound = ${dt_max}`, `optimal_iterations = 20`,
   `iteration_window = 0`, `growth_factor = 1.25`,
   `cutback_factor_at_failure = 0.5`, `timestep_limiting_function = pulse_sync`,
   `force_step_every_function_point = true`, `post_function_sync_dt = ${dt0}`,
   `timestep_limiting_postprocessor = dt_limit_pp`; postprocessor
   `dt_limit_pp` of type `FunctionValuePostprocessor` on `dt_limit`,
   `execute_on = 'INITIAL TIMESTEP_END'`.
8. Energy postprocessors: `external_work` gets `boundary = borehole`; the
   three `damping_work*` postprocessors and the `fdamp*` aux variables are
   removed; `full_input_energy = -external_work`; the `Receiver`, kinetic,
   elastic and `full_energy` blocks stay.
9. New measures: aux `r_damaged` (Lagrange, first) by `ParsedAux` with
   `coupled_variables = d`, `use_xyzt = true`,
   `expression = 'if(d >= ${d_th}, sqrt(x*x + y*y), 0)'`; aux `is_damaged`
   (constant monomial) with `expression = 'if(d >= ${d_th}, 1, 0)'`; a twin
   pair `r_damaged_090` and `fracture_extent_090` at the fixed threshold 0.9;
   postprocessors `r_dmg_max` (`NodalExtremeValue`, `value_type = max`),
   `fracture_extent` (`ParsedPostprocessor`, `max(r_dmg_max - ${rb}, 0)`),
   `damaged_volume` (`ElementIntegralVariablePostprocessor` of `is_damaged`),
   `crack_area` (`dissipated_energy_total / ${Gc_const}`),
   `p_wall_max` (`ElementExtremeValue` of `pulse_load_aux`, max),
   `borehole_area` (`AreaPostprocessor`, `boundary = borehole`, `execute_on = INITIAL`).
10. Outputs: Exodus `time_step_interval = 40`, `show = 'd pulse_load_aux'`,
    `execute_on = 'INITIAL TIMESTEP_END FINAL'`; CSV every step with the energy
    postprocessors and the new measures; `Checkpoint` every 100 steps, `num_files = 2`.
11. `fracture.i`: existing variables, bounds, kernels and materials;
    `initial_damage_aux` by one `ConstantAux` with `value = 0` and no block
    restriction; `Gc_var` and the Weibull comment removed;
    `cli_args = 'Gc_const=${Gc_const};l=${l};mesh_file=${mesh_file}'` in the
    `MultiApps` block.

**Interfaces.** Command-line overrides used by Phases 4 and 5: `mesh_file`,
`Np`, `eta`, `Gc_const`, `l`, `shape_alpha`, `shape_beta`, `zeta_pulse`,
`zeta_quench`, `t_free`, `Outputs/file_base`.

**Edge cases.**

- $N_p > 40$ $\to$ extend the `pulse_sync` and `dt_limit` lists.
- Nodes on the borehole edges at `zmin` and `zmax` carry both the pressure and the $u_z$ constraint; the Dirichlet condition wins for `disp_z`, as intended.

**Acceptance criteria.**

- [ ] `farms-opt -i elasticity.i mesh_file=../mesh/spe_cube_smoke.msh --check-input` passes.
- [ ] `farms-opt -i elasticity.i --mesh-only` lists the seven boundary names and block `acrylic`.
- [ ] Every header value equals the corresponding entry of `params_table.md`; `--show-input` prints `Gc_const` within 0.1 % of 2721 and the sub-app receives the same value.
- [ ] `grep -c FarmsNonReflectDashpotBC elasticity.i` prints 0, `grep -n "boundary = 'zmin zmax'" elasticity.i` hits the `disp_z` Dirichlet block, and `grep -n 'boundary = borehole' elasticity.i` hits both the pressure block and `external_work`.

**Depends on:** [Phase\ 1](#phase-1), [Phase\ 2](#phase-2)

## Phase 4: Smoke test on the coarse mesh {#phase-4}

**In one sentence.** Runs one pulse on the coarse mesh, undamped and damped,
and checks that the pressure, the energies, the sideset areas and the extent
metric behave as the formulation says.

**Implements [Data Flow](#data-flow) stage(s):** 9

**Files.** `code_3d/smoke/run_smoke.sh` (new) -- the two command lines.
`code_3d/smoke/check_smoke.py` (new) -- reads the CSVs and writes
`code_3d/smoke/SMOKE_REPORT.md`.

**Requirements.**

1. Run A (undamped):
   `mpiexec -n 8 farms-opt -i ../elasticity.i mesh_file=../../mesh/spe_cube_smoke.msh l=3e-3 Np=1 zeta_quench=0 Outputs/file_base=smokeA`
   (locally with the JIT shim on `PATH`, or on one Frontera node).
2. Run B (quench): as A with the default `zeta_quench=1e-5 Outputs/file_base=smokeB`.
3. `def check(csv_path: str, expected: dict) -> dict` computes
   (a) `borehole_area` against $2\pi r_b L_z$;
   (b) the peak of `p_wall_max` against \eqref{eq:peak} at $D = r_b$;
   (c) $\max_t|E_{tot}-E_{in}|/\max_t E_{in}$ for A, \eqref{eq:budget};
   (d) `external_work` nonzero in both runs;
   (e) `fracture_extent` and `crack_area` non-negative and non-decreasing;
   (f) $\mathcal{K}(T_w)/\max_t\mathcal{K}$ for A and B;
   (g) steps taken and wall time, to size [Phase\ 5](#phase-5).
4. `SMOKE_REPORT.md` lists each check with its threshold: (a) 0.5 %, (b) 2 %,
   (c) $\epsilon = 5$ %, (d) and (e) exact, (f) and (g) reported only.

**Interfaces.** The report's step count and wall time size [Phase\ 5](#phase-5).

**Edge cases.**

- The coarse mesh with $l = 3$ mm has $G_c = 13.6$ kJ/m$^2$ by \eqref{eq:strength} and a crack band wider than the borehole radius; this is a mechanics check, not physics.
- Check (c) above 5 % $\to$ rerun A with `Executioner/fixed_point_rel_tol=1e-8` before concluding a bookkeeping bug.
- Local JIT failure despite the shim $\to$ run on Frontera ([T1]).

**Acceptance criteria.**

- [ ] Checks (a) to (e) pass for runs A and B.
- [ ] $\mathcal{K}(T_w)/\max_t\mathcal{K}$ of run B is smaller than that of run A.
- [ ] Steps per pulse $\le 800$.

**Depends on:** [Phase\ 3](#phase-3)

## Phase 5: Production runs and comparison {#phase-5}

**In one sentence.** Submits each case to Frontera as one continuous pulse
train, continues it with `--recover` past each allocation limit, and plots the
extent after every completed pulse against the paper's four measurements.

**Implements [Data Flow](#data-flow) stage(s):** 10--11

**Files.** `code_3d/make_cases.py` (new) -- writes case directories.
`code_3d/case_template/submit_elasticity.sbatch`, `submit_elasticity_recover.sbatch`
(new, from `3d_solid/em_5J`) -- templates with placeholders. `code_3d/submit_all.sh`,
`code_3d/sync_results.sh` (copied from `3d_solid`, paths re-pointed).
`code_3d/post_fracture_extent.py` (new) -- the comparison.

**Requirements.**

1. `def write_case(root: Path, name: str, overrides: dict, nodes: int = 10, tasks: int = 500, hours: int = 48) -> Path`
   creates `root/name/` with `run.args` (one `key=value` per line) and the
   filled submit and recover scripts (job name `spe_<name>`, the
   `--split-mesh`/`--use-split` pair as in the existing scripts, the overrides
   appended to the `farms-opt` command).
2. Continuation: a job that reaches its wall-time limit is resubmitted with
   `submit_elasticity_recover.sbatch` (`--recover` from the latest checkpoint)
   until `end_time` or until the allocation is spent; the number of completed
   pulses is whatever the budget allows, the paper's twenty is not binding.
3. `def load_case(case_dir: Path, Tw: float) -> DataFrame` reads the CSV parts
   written by the successive runs, drops duplicated times keeping the last,
   and assigns the pulse index $k = \lfloor t/T_w\rfloor + 1$;
   `def extent_per_pulse(df, Tw) -> ndarray` takes `fracture_extent` at the
   last sample of each completed window. Figures: extent versus pulse number
   and versus cumulative stored energy $144\,k$ J with the points (5, 7.97),
   (11, 13.90), (15, 15.66), (20, 20.45) mm; $E_{in}$, $\mathcal{K}$, $\Psi_e$,
   $\Psi_d$ versus time; $A_f$ versus cumulative energy.
4. Run matrix written by `make_cases.py`; `base` runs first and for as many
   pulses as the budget allows, the others to at least five pulses:

| Case | $\eta$ | $G_c$ (J/m$^2$) | $t_d$ ($\mu$s) | $\zeta$ (s) | mesh |
| ---- | ------ | --------------- | -------------- | ----------- | ---- |
| `base` | 0.3 | 2721 | 6 | $2\times10^{-7}$ | h0p3 |
| `eta0p1` | 0.1 | 2721 | 6 | $2\times10^{-7}$ | h0p3 |
| `eta1` | 1.0 | 2721 | 6 | $2\times10^{-7}$ | h0p3 |
| `gc350` | 0.3 | 350 (literature, $\sigma_c = 24.7$ MPa) | 6 | $2\times10^{-7}$ | h0p3 |
| `td3` | 0.3 | 2721 | 3 | $2\times10^{-7}$ | h0p3 |
| `td12` | 0.3 | 2721 | 12 | $2\times10^{-7}$ | h0p3 |
| `zeta0` | 0.3 | 2721 | 6 | 0 | h0p3 |
| `eta1_weibull` | 1.0 | 2721, Weibull shape 10 | 6 | $2\times10^{-7}$ | h0p3 |
| `eta1_gauss` (2D first, Discovery) | 1.0 | 2721, Gaussian 12 % | 6 | quench | section h0p3 |

5. `sync_results.sh` pulls `*.csv` and `*.e` of every `cases/*` into the
   external-drive tree.

**Interfaces.** Overrides are the [Phase\ 3](#phase-3) header variables;
`post_fracture_extent.py` reads the CSV columns of its requirement 9.

**Edge cases.**

- Job dies or times out $\to$ `submit_elasticity_recover.sbatch`; no checkpoint found $\to$ resubmit the original script.
- The `gc350` override `Gc_const=350` leaves the header variable `ft` unused and MOOSE refuses unused parameters $\to$ its command line carries `--allow-unused`.
- A recover writes a new CSV file instead of appending $\to$ `load_case` merges the parts (requirement 3).

**Acceptance criteria.**

- [ ] On the smoke case, `--half-transient` followed by `--recover` reproduces the uninterrupted run's final CSV row within $10^{-8}$ relative.
- [ ] `submit_all.sh --dry-run` lists every `cases/*/submit_elasticity.sbatch`, each containing its `run.args` overrides.
- [ ] `analysis/fracture_extent.png` shows the four measured points and the `base` curve over the completed pulses.
- [ ] For `zeta0`, $|E_{tot}-E_{in}| \le 5$ % of $\max_t E_{in}$ over the first five pulses.

**Depends on:** [Phase\ 4](#phase-4)

## Testing Strategy

Phase 1: geometric invariants (areas, volume, counts, quality). Phase 2: the
script against the plan's tables and the production self-check. Phase 3:
`--check-input` and `--mesh-only`. Phase 4: energy conservation of the undamped
run, the pressure peak against \eqref{eq:peak}, irreversibility of the extent.
Phase 5: the half-transient recover test on the smoke case, then the
comparison plot. No closed-form fracture solution exists for this geometry;
the early-time hoop stress of Appendix A.8 is the only analytical anchor.

## Risks

| Risk | Severity | Mitigation |
| ------------------------------------------------------- | -------- | ----------------------------------- |
| Matching $f_t$ at $l = 0.6$ mm inflates $G_c$ to 2721 J/m$^2$, eight times the literature value, so cracks propagate too reluctantly and the calibrated $\eta$ is biased high. | high | report $G_c/350$ with every case; run `gc350`. |
| Twenty pulses cost about 15 000 steps on a 15 M mesh, fifteen times the existing production run. | high | checkpoint and `--recover` chain; stop when the budget is spent (the count is not binding); the `h0p4` mesh for sweeps. |
| The peak-pressure law was fitted for rock and water, not acrylic and brine. | high | treat $\eta$ as the calibration knob; bracket with the `peak_pressure` override. |
| Waves reflected by the free side faces and the clamped ends keep feeding damage after the pulse. | med | $\zeta$ sensitivity (`zeta0` case); 40 $\mu$s window sweep on the coarse mesh. |
| A `--recover` chain of a 17 M-element job may drop steps or split the CSV at each resubmission. | low | the half-transient recover test; `load_case` merges CSV parts. |

## Sources

\begingroup\footnotesize

| Id | What it is | Where it lives |
| ---- | ---------------------------------------- | ---------------------------------------------- |
| P1 | SPE-HFTC manuscript draft: Table 1 ($\rho$, $E$, $\nu$, $f_t$, compressive strength of cast acrylic), Table 2 (geometry, 144 J, 3--12 $\mu$s, 230--350 ns, 20 pulses at 80 s), Fig. 1 extents 7.97, 13.90, 15.66, 20.45 mm after 5, 11, 15, 20 pulses | `/Users/chunhuizhao/Downloads/` `SPE-HFTC-Manuscript-2026.pdf` |
| S1 | RACCOON spectral-split elasticity: stress, $\psi^+$, $\psi_e$ formulas; read 2026-10-21 | `src/materials/phasefield_smalldeform/` `SmallDeformationIsotropicElasticity.C`, `ComputeSmallDeformationStress.C` |
| S2 | Pulse function: distance to the gap segment in mm, peak-pressure law in bar and kJ, double exponential with $t \bmod T_w$ | `src/functions/pulsepower/` `ElkPulseLoadExperiment.C` |
| S4 | Nodal, boundary-restricted external work; the 3D inputs restrict it to the wrong surface (`'4'` with the pressure forces of `'3'`, and the reverse for damping) | `src/postprocessors/` `FarmsExternalWork.C`, `FirstStepElementIntegral` `VariablePostprocessor.C`; `3d_solid/em_5J/elasticity.i` lines 517--531 |
| S5 | Adaptive timestepper: growth cap, `force_step_every_function_point`, `post_function_sync_dt`, `timestep_limiting_postprocessor` | `src/timesteppers/` `FarmsIterationAdaptiveDT.C` |
| S6 | Production mesh: 2.64 M nodes, 16.37 M elements at 0.2 mm on 28 970 mm$^3$, giving 4.52 tetrahedra per $h^3$ | `pulsepower/cmame_revision/3d_mesh/` `cylinder_sample.msh` header and `.geo` |
| S7 | Existing production inputs and Frontera scripts, the template of this plan | `pulsepower/cmame_revision/` `3d_solid/em_5J/` |
| S8 | RACCOON phase-field kernels and functions: diffusion and source residuals, $c_0$ integration, degradation function, $\psi_d$, bounds | `src/kernels/phasefield/` `ADPFFDiffusion.C`, `ADPFFSource.C`; `src/materials/phasefield_function/*.C`; `src/bounds/` `VariableConstantIrreversibleBounds.C` |
| M1 | MOOSE dynamic kernels: the $\zeta\dot{\boldsymbol\sigma}$ Rayleigh term and $\rho\mathbf{a}$; kinetic energy $\tfrac12\rho|\mathbf{v}|^2$ | `raccoon/moose/modules/solid_mechanics/src/` `kernels/ADDynamicStressDivergenceTensors.C`, `kernels/InertialForce.C`, `auxkernels/KineticEnergyAux.C` |
| HM12 | PMMA parameters used in dynamic phase-field literature, $\rho = 1190$ kg/m$^3$, $E = 3.24$ GPa, $G_c = 350$ J/m$^2$ (Hofacker and Miehe 2012), as recorded in the repo's rate-dependent plan; the original paper was not re-read | `pulsepower/pf_puresolid_ratedependent/` `IMPLEMENTATION_PLAN.md`, section 3.1 |
| E1 | Exploration notes for the existing setup: call graph, conventions, gotchas | `pulsepower/spe_hftc/` `EXPLORE_3d_solid.md` |
| T1 | Local toolchain note: `farms-opt` links with the deferred-symbol workaround; the FParser JIT links the same way through the `mpicxx` shim in `tools/jit_shim` (AD parsed materials refuse `enable_jit = false`) | memory note `toolchain-macos-` `tahoe-link-failure`; `tools/jit_shim/mpicxx` |

\endgroup

## Appendix A: Derivations

### A.1 Equation \eqref{eq:pff}

Start from the regularised energy implemented by [S8] and [S1]:
$$\Pi[\mathbf{u},d] = \int_\Omega \Big[g(d)\,\psi^+(\boldsymbol\varepsilon) + \psi^-(\boldsymbol\varepsilon)\Big] dV
+ \int_\Omega \frac{G_c}{c_0}\Big(\frac{\alpha(d)}{l} + l\,|\nabla d|^2\Big) dV .$$
The sub-app assembles `ADPFFSource` $= (w, \partial\psi/\partial d)$ with
$\psi = \alpha G_c/(c_0 l) + g\,\psi^+$ and `ADPFFDiffusion`
$= (\nabla w,\ 2 G_c l/c_0\ \nabla d)$. Summing them is the first variation of
$\Pi$ with respect to $d$; integrating the gradient term by parts with a natural
boundary condition gives
$$\frac{G_c}{c_0}\Big(\frac{\alpha'(d)}{l} - 2l\,\Delta d\Big) + g'(d)\,\psi^+ = 0 .$$
For AT1, $\alpha(d) = d$, so $\alpha'(d) = 1$. `CrackGeometricFunction`
computes $c_0 = 4\int_0^1\sqrt{\alpha(s)}\,ds = 4\cdot\tfrac{2}{3} = \tfrac{8}{3}$
numerically. Irreversibility is imposed as the bound $d \ge d_{\text{old}}$
by the variational-inequality solver, and $d \le 1$ by `ConstantBounds`.

### A.2 Equation \eqref{eq:strength}

Uniaxial tension, homogeneous $d$, so $\Delta d = 0$ and $\psi^+ = \tfrac12 E\varepsilon^2$.
Damage starts when the source term first vanishes at $d = 0$:
$$g'(0)\,\tfrac12 E\varepsilon^2 + \frac{G_c}{c_0 l} = 0,\qquad g'(0) = -2(1-\eta_r)\approx -2 ,$$
hence $E\varepsilon^2 = G_c/(c_0 l)$ and
$\sigma_c = E\varepsilon = \sqrt{E G_c/(c_0 l)} = \sqrt{3EG_c/(8l)}$ with $c_0 = 8/3$.
Inverting, $G_c = c_0\,l\,\sigma_c^2/E = 8\,l\,\sigma_c^2/(3E)$. The plan sets
$\sigma_c = f_t = 69$ MPa (user decision, 2026-10-21), so $G_c$ follows from the
length scale ($E = 2.80$ GPa):

| $l$ | $G_c = 8lf_t^2/(3E)$ | ratio to 350 J/m$^2$ | $l_{ch} = EG_c/f_t^2 = 8l/3$ | $\sqrt{EG_c}$ |
| --- | -------------------- | -------------------- | ---------------------------- | ------------- |
| 0.3 mm | 1360 J/m$^2$ | 3.9 | 0.80 mm | 1.95 MPa m$^{1/2}$ |
| 0.5 mm | 2267 J/m$^2$ | 6.5 | 1.33 mm | 2.52 MPa m$^{1/2}$ |
| 0.6 mm | 2721 J/m$^2$ | 7.8 | 1.60 mm | 2.76 MPa m$^{1/2}$ |
| 0.8 mm | 3627 J/m$^2$ | 10.4 | 2.13 mm | 3.19 MPa m$^{1/2}$ |

The literature-toughness alternative (the `gc350` case) keeps $G_c$ physical
and lowers the strength instead:

| $G_c$ (J/m$^2$) | $l$ for $\sigma_c = 69$ MPa | $\sigma_c$ at $l = 0.3$ mm | at 0.6 mm | at 0.8 mm | $l_{ch}$ |
| --------------- | --------------------------- | -------------------------- | --------- | --------- | -------- |
| 250 | 55 $\mu$m | 29.6 MPa | 20.9 MPa | 18.1 MPa | 0.147 mm |
| 350 | 77 $\mu$m | 35.0 MPa | 24.7 MPa | 21.4 MPa | 0.206 mm |
| 500 | 110 $\mu$m | 41.8 MPa | 29.6 MPa | 25.6 MPa | 0.294 mm |

The existing granite case ($E = 50$ GPa, $G_c = 40$ J/m$^2$, $l = 0.2$ mm) has
$\sigma_c = 61$ MPa, which is why the same recipe cannot be copied: acrylic
is eighteen times softer, so matching 69 MPa needs a length scale eight times
smaller than the one used for granite.

### A.3 Equation \eqref{eq:shape}

With $f(\tau) = e^{-a\tau} - e^{-b\tau}$, $f'(\tau) = 0$ gives $a e^{-a\tau} = b e^{-b\tau}$,
so $\tau_p = \ln(b/a)/(b-a)$. The normalisation in \eqref{eq:pulse} makes
$p(t_r) = P_p$; the peak equals $P_p$ only when $\tau_p = t_r$, which is the
first condition. The second condition fixes the decay. Solutions for
$t_r = 0.3\ \mu$s (computed with `scipy.optimize.brentq` on 2026-10-21):

| $t_d$ ($\mu$s) | $a$ (1/s) | $b$ (1/s) | $p(40\ \mu\text{s})/P_p$ |
| -------------- | --------- | --------- | ------------------------ |
| 3 | $8.949\times10^{5}$ | $8.332\times10^{6}$ | $4\times10^{-16}$ |
| 6 | $4.103\times10^{5}$ | $1.153\times10^{7}$ | $9\times10^{-8}$ |
| 12 | $1.980\times10^{5}$ | $1.451\times10^{7}$ | $4\times10^{-4}$ |

The existing production parameters ($a = 4.658\times10^5$, $b = 4.661\times10^5$,
$t_r = 3\ \mu$s) do not satisfy the first condition: their peak is at
2.15 $\mu$s and reaches $1.065\,P_p$; the pulse decays to 10 % at 10.5 $\mu$s.

### A.4 Equation \eqref{eq:rayleigh}

The kernel [M1] adds $\zeta\,\dot{\boldsymbol\sigma}$ to the stress, i.e. a
damping force $\zeta\,\mathbf{K}\mathbf{v}$ for a linear structure with stiffness
$\mathbf{K}$. For a mode of stiffness $k_m$, mass $m_0$ and $\omega^2 = k_m/m_0$,
the damping constant is $c = \zeta k_m$ and the ratio is
$\xi = c/(2 m_0\omega) = \zeta k_m/(2 m_0 \omega) = \tfrac12\zeta\omega$.
A target of $\xi = 3$ % at 50 kHz gives $\zeta = 2\xi/\omega = 1.9\times10^{-7}$ s;
the plan uses $2\times10^{-7}$ s, i.e. 1.6 % at 25 kHz, 3.1 % at 50 kHz and
6.3 % at 100 kHz. The choice of 3 % is an assumption about acrylic, listed in
Open Questions.

### A.5 Mesh size estimates

From [S6], the production mesh has $N_e = 4.52\,V/h^3$. The fine zone is the
cube cross-section inside the disc $r \le r_f = 26$ mm (which touches all four
faces, so a crack reaching the face at $r = 24.4$ mm stays fine) minus the
borehole: $2118 - 54.1 = 2064$ mm$^2$, times $L_z = 49.14$ mm gives
$1.014\times10^{5}$ mm$^3$. The four corner regions hold $2.6\times10^{4}$ mm$^3$.

| Configuration | $h_f$ | fine zone | corners at 1 mm | total |
| ------------- | -------------- | --------- | --------------- | -------------------- |
| `h0p25` | 0.25 mm | 29.3 M | 0.12 M | 29.4 M (over budget) |
| `h0p3` | 0.30 mm | 17.0 M | 0.12 M | 17.1 M |
| `h0p4` | 0.40 mm | 7.2 M | 0.12 M | 7.3 M |
| `smoke` | 1.5 mm uniform | -- | -- | 0.17 M |

The meshes built with HXT on 2026-10-21 have 14.6 M (`h0p3`, fine zone
$r_f = 28$ mm so that the face centres stay fine) and 0.15 M (`smoke`)
tetrahedra, within the 20 % tolerance of these estimates.

### A.6 Peak pressures from \eqref{eq:peak}

$E_M = 0.144$ kJ, $B = 8000$, $\kappa = 0.35$, $m = 0.25$; $D = r_b = 4.15$ mm at
the gap. With the gap centred at $L_z/3 = 16.38$ mm, the near borehole end is at
$D = \sqrt{r_b^2 + (L_z/3 - L_g/2)^2} = 13.45$ mm and the far end at
$D = \sqrt{r_b^2 + (2L_z/3 - L_g/2)^2} = 29.47$ mm.

| $\eta$ | $P_p$ at the gap | near end | far end | with $B = 9000$, $m = 1$ (library defaults) |
| ------ | ---------------- | -------- | ------- | ------------------------------------------- |
| 1.0 | 284.4 MPa | 212.0 MPa | 174.3 MPa | 110.1 MPa |
| 0.5 | 223.2 MPa | 166.3 MPa | 136.7 MPa | 86.3 MPa |
| 0.3 | 186.6 MPa | 139.1 MPa | 114.3 MPa | 72.2 MPa |
| 0.1 | 127.1 MPa | 94.7 MPa | 77.8 MPa | 49.2 MPa |

For comparison the existing 5 J granite case has $P_p = 111$ MPa at its wall
($D = 1.6$ mm). Every row is above $\sigma_c = f_t = 69$ MPa, so damage starts
at the wall during the rise in all cases.

### A.7 Step budget per pulse

Growth from $\Delta t_0 = 10^{-8}$ s by 1.25 per step reaches $\Delta t_r = 3\times10^{-8}$ s
after 5 steps (0.08 $\mu$s); holding $\Delta t_r$ to 0.6 $\mu$s takes 17 steps;
growth to $\Delta t_{max} = 10^{-7}$ s takes 5 more (0.3 $\mu$s); the rest of the
60 $\mu$s window takes 591 steps. That is 618 steps per pulse without cutbacks,
about 740 with 20 % cutbacks and 15 000 for 20 pulses. The 2D section run of
2026-10-21 followed exactly this schedule (80 steps for the first 2 $\mu$s).
The existing production run (100 $\mu$s at $\Delta t_{max} = 10^{-7}$ s) takes
about 1000 steps.

### A.8 Early-time hoop stress

For a pressurised circular hole of radius $r_b$ in an infinite linear-elastic
medium (Lamé), $\sigma_{\theta\theta}(r) = p\,r_b^2/r^2$, tensile, equal to $p$
at the wall. Quasi-statically the hoop stress exceeds $\sigma_c$ out to
$r = r_b\sqrt{P_p/\sigma_c}$, i.e. 6.8 mm for $P_p = 187$ MPa and
$\sigma_c = 69$ MPa (5.6 mm at $\eta = 0.1$, 8.4 mm at $\eta = 1$; 11.4 mm for
the `gc350` strength of 24.7 MPa): an order-of-magnitude expectation for the
crushed-zone radius of the model, to be compared with the "intense, closely
spaced damage" zone the paper describes.

## Appendix B: Rationale

### Why damp at all, and why not a constant Rayleigh coefficient?

The experiment leaves 80 s between pulses, so each pulse meets a specimen at
rest. The model leaves 40 to 60 $\mu$s and resets nothing (user decision). Its
only energy sinks are the crack energy and the residual stiffness $\eta_r$, and
Newmark average acceleration conserves the rest. A 1D impedance estimate at
$\eta = 0.3$ ($P_p = 187$ MPa, $t_d = 6\ \mu$s, $\rho c_p = 2.3\times10^{6}$
Pa s/m, wall area $1.28\times10^{-3}$ m$^2$) puts the work done by one pulse at
30 to 40 J, of which the cracks of one pulse take about 1 J (300 mm$^2$ at
2721 J/m$^2$). Undamped, some 97 % of every pulse stays as vibration; after
twenty pulses the 0.15 kg cube would carry several hundred joules, a root-mean-
square velocity near 100 m/s and dynamic stresses of order 200 MPa everywhere,
and the model would disintegrate for a reason the experiment never sees.
Acrylic has a loss factor of a few per cent at these frequencies and the
holder absorbs more, so the real specimen is quiet long before the next pulse.

A constant stiffness-proportional coefficient cannot do this job cleanly,
because its damping ratio grows with frequency, \eqref{eq:rayleigh}. With
$\zeta = 2\times10^{-7}$ s the ring-down of the 20 to 100 kHz cube modes takes
100 $\mu$s or more, while the 0.3 $\mu$s rise, whose content sits near 1.7 MHz,
is critically damped ($\xi \approx 1$) and the damping stress during the rise
is $\zeta/t_r = 0.67$ of the elastic stress. Mass-proportional damping has the
opposite defect: it damps the low modes but not the stress waves that drive
cracks. The recommended alternative is a quench: $\zeta = 0$ for the first
25 $\mu$s of each window (the pulse, the primary crack growth, and the first
reflection, which returns to the borehole after 22 $\mu$s), then
$\zeta = 10^{-5}$ s for the remaining 35 $\mu$s. At that value every mode above
32 kHz is over-damped with a decay time $1/\zeta = 10\ \mu$s and the 20 kHz
modes decay in 13 $\mu$s, so the window ends with under 1 % of the vibration
energy and the damaged cube relaxed to static equilibrium, which is the state
the experiment is in before the next pulse. It needs no new code: `zeta` is a
material-property name, so a `GenericFunctionMaterial` driven by a
`PiecewiseConstant` time function supplies the switch. The price is a 60 $\mu$s
window (about 620 steps per pulse instead of 420 before cutbacks) and an
energy budget whose gap equals the quenched energy.

### How $\zeta$ enters the equations

*Continuum.* The kernel solves \eqref{eq:momentum} with the stress replaced by
$\boldsymbol\sigma + \boldsymbol\sigma_v$, where $\boldsymbol\sigma_v = \zeta\,\dot{\boldsymbol\sigma}$.
Because $\boldsymbol\sigma = \mathbb{C}(d)\!:\!\boldsymbol\varepsilon$ with the degraded
tensor of \eqref{eq:stress}, the damping stress is $\zeta\,\mathbb{C}(d)\!:\!\dot{\boldsymbol\varepsilon}$:
a Kelvin--Voigt solid whose viscosity tensor is $\zeta$ times its (degraded)
stiffness. $\zeta$ therefore has the unit of time, and the material is
"as viscous as it is stiff" over a time scale $\zeta$. The dissipated power
density is $\zeta\,\dot{\boldsymbol\varepsilon}\!:\!\mathbb{C}\!:\!\dot{\boldsymbol\varepsilon} \ge 0$,
and its time integral is the gap in \eqref{eq:budget}.

*Discrete.* After finite-element assembly the semi-discrete system is
$\mathbf{M}\ddot{\mathbf{u}} + \mathbf{C}_d\dot{\mathbf{u}} + \mathbf{K}\mathbf{u} = \mathbf{f}(t)$
with $\mathbf{C}_d = \zeta\,\mathbf{K}$: stiffness-proportional Rayleigh damping
(the mass-proportional part, `eta` of `ADInertialForce`, is not used). For a
mode of natural frequency $\omega$ the modal equation is
$\ddot q + \zeta\omega^2\dot q + \omega^2 q = 0$, so the damping ratio is
$\xi = \tfrac12\zeta\omega$ (equation \eqref{eq:rayleigh}) and the loss factor
is $\tan\delta = \zeta\omega$: damping grows linearly with frequency.

*Kernel.* `ADDynamicStressDivergenceTensors` approximates the stress rate by a
backward difference, $\dot{\boldsymbol\sigma} \approx (\boldsymbol\sigma_{n+1}-\boldsymbol\sigma_n)/\Delta t$,
and with the HHT parameter at zero assembles
$$\int_\Omega \nabla w \cdot \Big[\big(1+\tfrac{\zeta}{\Delta t}\big)\boldsymbol\sigma_{n+1} - \tfrac{\zeta}{\Delta t}\,\boldsymbol\sigma_n\Big]\,dV ,$$
read from [M1]. `zeta` is a material property evaluated at $t_{n+1}$, which is
why a time function (`zeta_fn` through `GenericFunctionMaterial`) can switch it
within a window. The stress at $t_n$ is the stored old state, so the term costs
nothing extra.

*Numbers.* Two regimes follow from $\xi = \tfrac12\zeta\omega$:

- Under-damped modes ($\zeta\omega < 2$) decay with time constant
  $1/(\xi\omega) = 2/(\zeta\omega^2)$.
- Over-damped modes ($\zeta\omega > 2$) do not oscillate; their slow root is
  $-1/\zeta$, so they decay with time constant $\zeta$ whatever their frequency.

With $\zeta = 10^{-5}$ s every mode above $1/(\pi\zeta) = 32$ kHz is over-damped
and decays in 10 $\mu$s, and the 20 kHz cube modes decay in 13 $\mu$s, so a
35 $\mu$s quench removes over 99 % of the vibration energy. Because the damping
force vanishes with the velocity, the quench does not alter the static
equilibrium the damaged cube settles into. With $\zeta = 2\times10^{-7}$ s held
during the pulse, the strain rate of the 0.3 $\mu$s rise gives a damping stress
of $\zeta/t_r = 0.67$ of the elastic stress and $\tan\delta = 2$ at 1.7 MHz,
forty times the loss factor of real acrylic there, which is why the schedule
sets $\zeta = 0$ while the pulse acts.

### Why back-calculate $G_c$ from $f_t$ instead of using the literature toughness?

Decision taken by the user on 2026-10-21 after the first draft of this plan.
In AT1 the three quantities $\sigma_c$, $G_c$ and $l$ are tied by
\eqref{eq:strength}, and the mesh fixes $l \ge 2h_f$. At the affordable
$l = 0.6$ mm one of the two material numbers must be wrong: either the
strength drops to 24.7 MPa with $G_c = 350$ J/m$^2$, or the toughness rises to
2721 J/m$^2$ with $\sigma_c = 69$ MPa. The two errors act on different parts of
the physics. The strength controls where damage nucleates and how large the
crushed zone around the borehole becomes (Appendix A.8: 6.8 mm radius at
69 MPa versus 11.4 mm at 24.7 MPa for the baseline pressure). The toughness
controls how far a crack can run for a given energy input, and therefore the
extent per pulse that the paper measures. The user chose to keep the measured
strength, so that nucleation and the crushed zone are physical, and to let the
energy-efficiency calibration absorb the inflated toughness: a larger $\eta$
delivers more energy and lengthens the cracks. The choice was confirmed on
2026-10-21 after weighing a finer length scale: with $h_f = l/2$ the element
count grows as $l^{-3}$ and the step count as $l^{-1}$, so $l = 0.3$ mm would
cost about sixteen times the planned run for an inflation of 3.9 instead of
7.8, and $l = 77\ \mu$m is out of reach. Consequences to keep in view:

- The calibrated $\eta$ is biased high relative to the true energy share;
  report it together with $G_c/350 = 7.8$.
- The crack-area estimate $A_f = \Psi_d/G_c$ of \eqref{eq:area} divides by the
  inflated $G_c$; multiply by 7.8 to compare with a 350 J/m$^2$ material.
- The inflation is linear in $l$: $l = 0.3$ mm would reduce it to 3.9 at eight
  times the element count (Appendix A.5), which is not affordable in 3D.
- An alternative would have been the PF-CZM degradation function, which
  matches $f_t$ with the literature $G_c$ for any resolved $l$; the user chose
  to stay with AT1 (2026-10-21), so it is no longer in the plan.
- The efficiency knob cannot fully offset the inflation: $P_p \propto \eta^{0.35}$,
  so raising $\eta$ from 0.3 to its ceiling of 1 raises the pressure 1.5 times
  and the crack driving force about 2.3 times, against the 7.8 demanded. If
  `base` underpredicts the measured extents by more than that, the toughness
  is the limiting factor. The fracture energy itself is not: all cracks after
  20 pulses cost of order 10 to 20 J against 2880 J stored.
- The `gc350` case of the run matrix keeps the alternative, so the two choices
  bracket the calibrated extent.
