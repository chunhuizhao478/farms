# `pulsepower/cmame_revision/3d_solid` — 3D pure-solid phase-field pulse-power setup

Exploration written 2026-10-21 as groundwork for adapting this setup to the
SPE-HFTC-2026 acrylic-cube experiment (see `PLAN.md` in this folder).

## Overview

`3d_solid` is a pair of MOOSE input files (`elasticity.i` + `fracture.i`) run
with the `farms-opt` application. It simulates a hollow granite cylinder loaded
on its borehole wall by a train of electro-hydraulic pressure pulses, with
dynamic (Newmark) elasticity in the main app and an AT1 phase-field fracture
sub-app coupled by fixed-point iteration. Two cases exist, `em_2d5J` (2.5 J per
pulse) and `em_5J` (5 J per pulse); they differ only in the `EM` parameter of
the pulse function. There is no C++ in this folder; all physics lives in
`src/` (farms), `raccoon/` (RACCOON phase-field objects compiled into farms) and
`raccoon/moose/modules/solid_mechanics` (MOOSE).

## Scope

- Target: `pulsepower/cmame_revision/3d_solid/{em_2d5J,em_5J}/{elasticity.i,fracture.i}`,
  `submit_*.sbatch`, `submit_all.sh`, `sync_results.sh`, and the mesh
  `pulsepower/cmame_revision/3d_mesh/cylinder_sample.{geo,msh}`.
- Entry points: `elasticity.i` (main app, run with `farms-opt -i elasticity.i`);
  `fracture.i` is launched as a `TransientMultiApp` sub-app.
- Size: 2 input files (~620 + ~210 lines) per case; 14 custom farms/RACCOON
  objects referenced (listed in the walkthrough), plus MOOSE solid-mechanics
  kernels/aux kernels.
- Mesh: gmsh 4.1 ASCII, 2,638,532 nodes, 16,365,590 elements (tets + surface
  triangles), uniform characteristic length 0.2 mm; cylinder radius 12.5 mm,
  borehole radius 1.6 mm, height 60 mm. Physical groups: volume 1 "Pierced
  Cylinder"; surfaces 2 Bottom, 3 Inner (borehole), 4 Outer, 5 Upper.

## Architecture

- Key objects
  - `ElkPulseLoadExperiment` (farms Function): borehole pressure $p(t,\mathbf{x})$.
  - `Pressure` action (MOOSE): applies $-p\,\mathbf{n}$ on boundary 3, saves
    nodal reaction forces to `fx fy fz`.
  - `ADDynamicStressDivergenceTensors` + `ADInertialForce` (MOOSE): dynamic
    momentum balance; `NewmarkBeta` integrator with `NewmarkAccelAux/VelAux`.
  - `SmallDeformationIsotropicElasticity` (RACCOON, `decomposition = SPECTRAL`)
    + `ComputeSmallDeformationStress` + `ADComputeSmallStrain`: degraded stress
    and the active strain energy `psie_active` that drives damage.
  - `PowerDegradationFunction` $g(d)=(1-d)^2(1-\eta_r)+\eta_r$, $\eta_r=10^{-6}$.
  - `FarmsNonReflectDashpotBC` (farms): Lysmer dashpots on boundaries 2, 4, 5.
  - `FarmsIterationAdaptiveDT` (farms): iteration-adaptive $\Delta t$ with a
    hard cap `max_time_step_bound`.
  - Sub-app: `ADPFFDiffusion` + `ADPFFSource` (RACCOON kernels), `CrackGeometricFunction`
    (AT1, $\alpha=d$, computes $c_0$), `ADDerivativeParsedMaterial` `psi`,
    `CrackDissipatedEnergyDensity` (farms), `VariableConstantIrreversibleBounds` +
    `ConstantBounds` (irreversibility via PETSc VI solver `vinewtonrsls`).
  - Energy bookkeeping: `FarmsExternalWork` (nodal, boundary-restricted),
    `ADKineticEnergyAux`, `ADElementIntegralMaterialProperty` (psie),
    `FirstStepElementIntegralVariablePostprocessor`, `Receiver` + `ParsedPostprocessor`.
- Data flow per time step:
  mesh + $p(t,\mathbf{x})$ $\to$ main app solves $\mathbf{u}$ with the current
  $d$ $\to$ `psie_active`, `mesh_size` copied to the sub-app $\to$ sub-app
  solves $d$ (bounded) $\to$ `d`, `dissipated_energy_density` copied back $\to$
  fixed-point loop until `fixed_point_rel_tol = 1e-6` $\to$ postprocessors
  (energies) $\to$ CSV every 40 steps, Exodus (`d`, `pulse_load_aux`) every 40
  steps, checkpoint every 100 steps.
- External dependencies: MOOSE solid_mechanics module (kernels, Newmark aux,
  Pressure action, KineticEnergyAux), PETSc (hypre BoomerAMG + GMRES,
  `vinewtonrsls`), gmsh (mesh), SLURM on Frontera (`--split-mesh 500`
  pre-partitioning then `--use-split`).

## Call Graph (critical path, one time step)

```
Transient executioner (elasticity.i)
  └→ FarmsIterationAdaptiveDT::computeDT()            dt from iteration count, capped at 1e-7 s
  └→ fixed-point loop (max 10 its)
       ├→ main app Newton solve (SMP, hypre/GMRES, vinewtonrsls)
       │    ├→ ADDynamicStressDivergenceTensors (x,y,z)   ∇w·σ  (+ HHT/Rayleigh terms, both 0 here)
       │    │    └→ ComputeSmallDeformationStress
       │    │         └→ SmallDeformationIsotropicElasticity::computeStressSpectralDecomposition
       │    │              σ = g(d) σ⁺ + σ⁻ ; psie_active = ψ⁺(ε) ; psie = g ψ⁺ + ψ⁻
       │    ├→ ADInertialForce (x,y,z)                   ρ a  (Newmark-β)
       │    ├→ Pressure BC on boundary 3                  −p(t,x) n, p from ElkPulseLoadExperiment::value
       │    └→ FarmsNonReflectDashpotBC on 2,4,5          ρ(c_p v_n n + c_s v_t)
       ├→ MultiAppCopyTransfer psie_active, mesh_size → fracture sub-app
       ├→ sub-app Newton solve (vinewtonrsls with Bounds)
       │    ├→ ADPFFDiffusion     2 G_c l / c_0 ∇w·∇d
       │    ├→ ADPFFSource        w ∂ψ/∂d,  ψ = α(d) G_c/(c_0 l) + g(d) psie_active
       │    └→ CrackDissipatedEnergyDensity → dissipated_energy_density (aux + PP integral)
       └→ MultiAppCopyTransfer d, dissipated_energy_density → main app
  └→ AuxKernels: NewmarkAccelAux/VelAux, FunctionAux(pulse_load_aux), ADKineticEnergyAux
  └→ Postprocessors: external_work, damping_work*, full_input_energy, energies, full_energy
```

## Detailed Walkthrough

### `elasticity.i` header parameters
- $E = 50$ GPa, $\nu = 0.373$, $G_c = 40$ J/m², $\rho = 2600$ kg/m³, $l = 0.2$ mm
  (comment: "l = 5e-4, ft = 61 MPa" is stale; with $l=0.2$ mm the AT1
  strength is $\sigma_c=\sqrt{3EG_c/(8l)} = 61$ MPa).
- Newmark $\beta=0.25$, $\gamma=0.5$ (average acceleration, no numerical
  damping), `hht_alpha = 0`.
- $c_p=\sqrt{(K+4G/3)/\rho}$, $c_s=\sqrt{G/\rho}$ feed the dashpots.

### `ElkPulseLoadExperiment::value(t, p)` (`src/functions/pulsepower/ElkPulseLoadExperiment.C`)
- Purpose: borehole-wall pressure as a function of time and position.
- Inputs: `shape_param_alpha` $a$, `shape_param_beta` $b$, `rise_time` $t_r$,
  `single_pulse_duration` $T_w$, `EM` (energy, **kJ**), `gap` $L_g$ (m),
  `convert_efficiency` $\eta$, `fitting_param_alpha` $\kappa$,
  `fitting_param_exponent` $m$ (default 1), `base_factor` $B$ (default 9000),
  `discharge_center` $\mathbf{x}_c$, `number_of_pulses` $N_p$, optional
  `peak_pressure` override and `minimum_applied_pressure` floor.
- Algorithm:
  1. $D$ (**mm**) = 1000 × distance from the point to the discharge *segment*
     (the vertical segment of length $L_g$ centred at $\mathbf{x}_c$): radial
     distance inside the gap's z-range, 3D distance to the nearest gap end outside it.
  2. $P_p = 10^5\,B\,(\eta\,E_M)^{\kappa}/D^{m}$ Pa (the $10^5$ converts bar to Pa).
     If `peak_pressure > 0` it replaces $P_p$ (then no spatial variation).
  3. $\tau = t \bmod T_w$; for $t \le N_p T_w$:
     $p = P_p\,[e^{-a\tau}-e^{-b\tau}]/[e^{-a t_r}-e^{-b t_r}]$, else 0.
  4. Optional floor `max(p, minimum_applied_pressure)`.
- Math notes: the double exponential peaks at $t_p=\ln(b/a)/(b-a)$, **not** at
  $t_r$; $t_r$ only sets where $p = P_p$. With the production values
  $a=4.658\times10^5$, $b=4.661\times10^5$, $t_r=3$ µs: $t_p = 2.15$ µs,
  $p(t_p) = 1.065\,P_p$, and $p$ drops to 10 % of peak at 10.5 µs. Pulses are
  back-to-back: `number_of_pulses = 10` with $T_w = 10$ µs gives a 100 µs train.
- Production peak at the borehole wall ($D = 1.6$ mm, $E_M = 0.005$ kJ, $\eta=1$,
  $B=8000$, $\kappa=0.35$, $m=0.25$): $P_p = 111$ MPa.

### `Pressure` action (MOOSE) on boundary 3
Applies traction $-p\,\mathbf{n}$ with $p$ = function value at the quadrature
point (so the spatial decay along the borehole is honoured). `save_in_disp_*`
accumulate the nodal residual contributions into `fx fy fz`; those are nonzero
only at nodes of boundary 3.

### `FarmsNonReflectDashpotBC` (`src/bcs/core/FarmsNonReflectDashpotBC.C`)
- Residual $\int w\,\rho\,(c_p v_n n_i + c_s (v_i - v_n n_i))\,dA$ with the
  velocity reconstructed from the Newmark update inside the BC (uses old
  velocity/acceleration aux variables). Applied on Outer (4), Upper (5), Bottom
  (2) in the production input: the finite cylinder is treated as embedded in an
  absorbing far field, not as a free specimen.

### `SmallDeformationIsotropicElasticity` (SPECTRAL)
$\lambda = K - 2G/3$, $\boldsymbol\varepsilon^+$ = positive spectral part,
$\sigma^+ = \lambda\langle\mathrm{tr}\,\varepsilon\rangle_+ \mathbf{I} + 2G\varepsilon^+$,
$\sigma = g\,\sigma^+ + (\sigma_{\text{intact}}-\sigma^+)$,
$\psi^+ = \tfrac12\lambda\langle\mathrm{tr}\,\varepsilon\rangle_+^2 + G\,\varepsilon^+\!:\!\varepsilon^+$,
`psie` $= g\psi^+ + \psi^-$, $\partial\psi_e/\partial d = g'(d)\psi^+$.
Only tension-like strain drives damage; the crushed zone around the borehole
therefore forms through the hoop (tensile) stress, not through the compressive
radial stress.

### Sub-app `fracture.i`
- Variable $d$ (Lagrange, first order) bounded in $[\,$`initial_damage_aux`$,1]$
  by PETSc VI; `initial_damage_aux` is 0 on both blocks in `3d_solid` (the
  `SubdomainBoundingBoxGenerator` block 2 is a vestigial "weak plane" slot; in
  `3d_lab_experiment` it is set to 0.9).
- AT1: $\alpha(d)=d$, `CrackGeometricFunction` integrates $c_0 = 4\int_0^1\sqrt{\alpha}\,ds = 8/3$
  numerically; $\xi=\alpha'(0)=1$.
- $\psi = \alpha G_c/(c_0 l) + g(d)\,$`psie_active`; kernels assemble
  $(\nabla w, 2G_c l/c_0 \nabla d) + (w, \partial\psi/\partial d) = 0$, i.e.
  the strong form $\tfrac{G_c}{c_0}\left(\tfrac{1}{l} - 2l\Delta d\right) + g'(d)\psi^+ = 0$ where active.
- `dissipated_energy_density` $= \tfrac{G_c}{c_0}\left(\tfrac{\alpha}{l} + l|\nabla d|^2\right)$,
  integrated in the sub-app (`dissipated_energy_dynamic`) and sent to the main
  app through `MultiAppPostprocessorTransfer`; the value at the first step is
  frozen by `FirstStepElementIntegralVariablePostprocessor` and subtracted
  (`dissipated_energy_total = dynamic − first_step`) so that pre-seeded damage
  does not count as dissipation.
- `Gc_var` aux + a commented Weibull `RandomIC` show how a heterogeneous
  toughness field was tried; it is off in production.

### Time stepping (`FarmsIterationAdaptiveDT`)
- `dt = 1e-8` initial, grows by 1.25 when Newton iterations < 20, cut by 0.5 on
  failure, capped by `max_time_step_bound = 1e-7` (applied only in the growth
  branch). Also supports `time_t/time_dt`, `timestep_limiting_function` with
  `force_step_every_function_point` (lands exactly on the points of a
  `PiecewiseLinear`) and `post_function_sync_dt` (dt right after such a point);
  these are unused in `3d_solid` but are exactly what a multi-pulse train with
  quiet windows needs.
- `end_time = 1e-4` s = 10 pulses × 10 µs.

### Energy bookkeeping (main app)
- `external_work` = Σ over boundary nodes of $\mathbf{f}\cdot\dot{\mathbf{u}}\,\Delta t$,
  accumulated (`FarmsExternalWork`, transient mode).
- `full_input_energy = −external_work − damping_work − damping_work_top − damping_work_bottom`
  (sign: `save_in` stores residual contributions, i.e. minus the applied force).
- `full_energy = kinetic + elastic(psie) + dissipated_total`; the energy budget
  check is `full_energy ≈ full_input_energy`.

### Job scripts
- `submit_elasticity.sbatch`: 10 nodes / 500 MPI tasks, 48 h, first
  `ibrun -n 1 farms-opt ... --split-mesh 500 --split-file foo.cpr`, then
  `ibrun farms-opt ... --use-split`. `submit_elasticity_recover.sbatch` adds
  `--recover` (checkpoint every 100 steps, 2 files kept).
- `submit_all.sh` finds every `submit_elasticity.sbatch` under the folder;
  `sync_results.sh` rsyncs `*.e`/`*.csv` from Frontera to an external drive.

## Conventions

- Units: SI in the inputs (m, s, Pa, kg), **except** inside
  `ElkPulseLoadExperiment` where `EM` is kJ and $D$ is mm (bar-based empirical law).
- Boundaries are referenced by numeric gmsh physical tags, with comments
  stating the intent; names ("Inner", "Outer", ...) exist in the mesh and would
  be safer.
- Parameters are duplicated between `elasticity.i` and `fracture.i` by hand
  (mesh file, box generators); `Gc_const` and `l` are passed by `cli_args`.
- Case folders differ by a single line; `submit_all.sh` discovers them.
- AD objects (`AD*`) throughout; the dashpot BC is non-AD (hand Jacobian).

## Gotchas

1. **Energy postprocessor boundaries are swapped in the 3D inputs.**
   `external_work` is restricted to `boundary = '4'` (Outer) but its forces
   `fx fy fz` are nonzero only on boundary 3 (the pressure surface);
   `damping_work` is restricted to `boundary = '3'` but `fdamp*` live on
   boundary 4. `FarmsExternalWork` is a `NodalPostprocessor` and only visits
   nodes of its boundary, so both terms evaluate to zero; `full_input_energy`
   in these runs contains only the top/bottom dashpot work. The 2D input
   (`2d_solid/no_weak_point`) has them right (`'3'` for external, `'1'` for
   damping). Same defect in `3d_lab_experiment/phase_field_pure_solid`.
   Confirm on a production CSV (`full_input_energy` ≪ `full_energy`).
2. `rise_time` is a normalisation point, not the time of the peak (see above);
   the real peak is 6.5 % above `Pp` for the production shape parameters.
3. Pulses are contiguous (`t mod T_w`); there is no quiet interval and the
   function does not reset the mechanical state between pulses.
4. $D$ is measured to the discharge **segment**, so every point of the borehole
   wall inside the gap sees the same $P_p$; outside the gap $P_p$ decays only as
   $D^{-0.25}$ (with `fitting_param_exponent = 0.25`), so the whole borehole is
   strongly loaded.
5. `mesh_size` (`MeshSize` aux) assumes equilateral tets:
   $h = (6\sqrt2 V_e)^{1/3}$; it is only transferred to the sub-app and not used
   by any production material.
6. Production uses $h = l = 0.2$ mm (one element per regularisation length);
   the usual recommendation is $h \le l/2$.
7. The corner node fixed by `ExtraNodesetGenerator` at the node closest to
   (0.1, 0.1, 0) pins all three displacement components at one point of the
   bottom rim; with absorbing boundaries this is redundant but harmless.
8. `fixed_point_max_its = 10` with `accept_on_max_fixed_point_iteration = false`:
   a step that does not converge the $u$–$d$ loop is rejected and $\Delta t$ cut.
9. Sub-app `psi` is an `ADDerivativeParsedMaterial` (JIT-compiled). On this Mac
   the JIT link fails (see memory `toolchain-macos-tahoe-link-failure`); local
   runs need `enable_jit = false` on that material.
10. The gmsh `.geo` hard-codes OCC surface ids (`{4},{5},{6},{7}`); they change
    if the geometry script changes, so a new mesh script should identify faces
    geometrically.
11. `ErrorPsiMeasure` (aux `eng_err`) requests the *old* `psie_active` and thus
    makes it stateful; it is diagnostic only.
12. Exodus output shows only `d` and `pulse_load_aux` (`show = ...`), so stress
    and displacement are not available for post-processing from production runs.

## Open Questions

- Origin of the peak-pressure law constants ($B = 8000$, $\kappa=0.35$,
  $m=0.25$; defaults 9000 and 1): not documented in the repo (no script in
  `pulsepower/pulseloading` produces them). Needs the user's reference.
- Whether the production energy CSVs indeed show `full_input_energy ≈ 0`
  (confirms gotcha 1); the local copies of the CSVs contain only the header row.
- Wall time per step on Frontera for the 16.4 M-element mesh (needed to size
  longer pulse trains); the `.o` logs are on Frontera, not in the repo.
