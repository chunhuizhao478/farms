# Implementation Plan: Strain-Based Porosity Update Option in `ElkPorousFlowDamagedPorosity`

## Overview
Add a selectable porosity-update law to the damaged-porosity material so the existing
damage-driven (bounded-maximum) formula can be swapped for the strain-based update of
Liu et al. (2024, CMAME 429:117165), eq. (40): `phi = phi_0 + eps_1`, where `eps_1` is the
maximum (most-tensile) principal strain. Both laws share the same lower/upper porosity
clamps and the same output property name, so a study can switch between
"bounded-maximum porosity" and "strain-based porosity" by changing one input parameter and
compare the effect with no other changes.

## Decisions locked for this pass
- **Reversibility:** the strain law is **reversible / paper-literal** — `phi = phi_0 + eps_1`
  evaluated instantaneously each step (porosity falls if strain relaxes). No max-over-history
  state. (The monotone variant is noted only as a future extension in Risk Assessment.)
- **Scope:** implement the **non-AD** material now (Phase 1) with its tests (Phase 3) and docs
  (Phase 4). **Phase 2 (AD twin) is deferred** — do not implement it in this pass; it is kept
  in the plan only so the eventual mirror stays consistent.

## Background: the math

The material currently writes the **damage-based** porosity (paper eq. 36, [38]):

```
g(d)     = (1 - d)^2                          # AT1/AT2 degradation, d = damage in [0,1]
phi_raw  = phi_0 + (1 - phi_0) * (1 - g(d))
phi      = clamp(phi_raw, lower_bound, upper_bound)
```

As `d -> 1` this saturates at `phi = 1`, capped by `porosity_upper_bound` — this is the
"bounded maximum porosity" baseline the user wants to compare against.

The **strain-based** update proposed by the paper is derived as:

- eq. (39): a line crack of aperture `omega` inside a quad element of edge `h_e` raises the
  cell-average porosity by `omega/h_e`, i.e. `phi_1 = phi_m + omega/h_e`.
- eq. (32): the aperture is `omega = h_e * eps_1`, with `eps_1` the maximum principal strain.
- Substituting gives **eq. (40):** `phi_1(eps) = phi_m + eps_1`.

So the new branch computes:

```
eps     = mechanical_strain (rank-two kinematic strain tensor)
eps_1   = max principal strain = largest eigenvalue of eps
phi_raw = phi_0 + eps_1                        # Liu et al. (2024) eq. (40)
phi     = clamp(phi_raw, lower_bound, upper_bound)
```

Notes that drive the design:
- The paper's phase-field `v` is intactness (`v=1` intact); this codebase uses `d = 1 - v`
  (damage, `d=0` intact). The existing `g(d) = (1-d)^2` already encodes that mapping; the
  strain branch does not use `d` at all.
- `eps_1` is the **largest** eigenvalue. `RankTwoTensor::symmetricEigenvalues` returns
  eigenvalues in **ascending** order (LAPACK `dsyev`), so `eps_1 = eigvals.back()`. This
  matches the paper's ordering `eps_1 >= eps_2 >= eps_3` (eq. 29).
- In 2D plane strain `eps_zz = 0`, so `eps_1 = max(eps_xx_principal, eps_yy_principal, 0) >= 0`.
  Under pure compression the strain branch therefore yields `phi = phi_0` (no spurious
  porosity loss) — a desirable property, not a bug.

## Constraints

- **Output property name is fixed.** The branch must keep declaring
  `PorousFlow_porosity_qp_damaged` (qp) / `PorousFlow_porosity_nodal_damaged` (nodal). All
  downstream consumers (`ElkPorousFlowDamagedBiotModulus` with `use_damaged_porosity=true`,
  the `porosity_aux` `MaterialRealAux`, kinetic-energy `ParsedAux` blocks) reference these
  names and must not change.
- **Existing default behavior must not change.** Inputs that omit the new parameter must
  produce bit-identical results to today (damage model). The new parameter defaults to
  `damage`.
- **Strain source.** Bind the kinematic strain via the `mechanical_strain` rank-two material
  property declared by `ComputeSmallStrain` (the same property the normal-strain permeability
  branch in `NDSmallDeformationIsotropicElasticity` binds; see its comment at
  `src/materials/phasefield_smalldeform_nonAD/NDSmallDeformationIsotropicElasticity.C:299-315`).
  The property name must be overridable via a `MaterialPropertyName` parameter.
- **qp-only for the strain branch.** `mechanical_strain` exists only at quadrature points.
  The strain branch must `mooseError` if the material is being evaluated nodally
  (`isNodal()`), with a clear message. (Production inputs use only the qp instance.)
- **Convention constraints.** Follow the existing file style: `registerMooseObject("farmsApp", ...)`,
  `std::clamp`, MOOSE `MooseEnum` for the model switch (mirror the
  `permeability_model = none|exponential|darcy_poiseuille|normal_strain` style already used in
  this codebase). Keep the non-AD and AD twins in sync.
- **No new external dependencies.** Eigenvalues come from the existing
  `RankTwoTensor::symmetricEigenvalues`.

---

## Phase 1: Strain-based option in the non-AD material (primary deliverable)

### Goal
After this phase, `[porosity_damaged] type = ElkPorousFlowDamagedPorosity` accepts
`porosity_update_model = damage|strain`; `strain` evaluates `phi = clamp(phi_0 + eps_1)` from
the maximum principal `mechanical_strain`, while `damage` (default) is unchanged.

### Files to Modify
- `include/materials/porousflowmatprops/ElkPorousFlowDamagedPorosity.h` — add the model enum,
  the strain-property pointer member, and `#include "RankTwoTensor.h"`.
- `src/materials/porousflowmatprops/ElkPorousFlowDamagedPorosity.C` — add the two new params,
  parse/guard them in the constructor, and branch in `computeQpProperties()`.

### Detailed Requirements

1. **Header (`ElkPorousFlowDamagedPorosity.h`)**
   - Add `#include "RankTwoTensor.h"` after `#include "Material.h"`.
   - Inside the class, add a scoped enum and members:
     ```cpp
     /// Which porosity-update law to evaluate
     enum class PorosityUpdateModel { DAMAGE, STRAIN };
     const PorosityUpdateModel _porosity_update_model;

     /// Kinematic strain tensor (bound only when _porosity_update_model == STRAIN; else nullptr).
     /// eps_1 = max principal value of this tensor drives the eq. (40) update.
     const MaterialProperty<RankTwoTensor> * _mechanical_strain;
     ```
   - Keep all existing members (`_damage`, `_initial_porosity`, `_porosity_lower_bound`,
     `_porosity_upper_bound`, `_porosity_damaged`) unchanged.
   - Update the class doxygen comment to document both laws.

2. **`validParams()` additions** (after the existing `porosity_upper_bound` param):
   ```cpp
   params.addParam<MooseEnum>(
       "porosity_update_model",
       MooseEnum("damage strain", "damage"),
       "Porosity update law. 'damage' (default): phi = phi_0 + (1-phi_0)[1-(1-d)^2], the "
       "existing damage-driven bounded-maximum-porosity model. 'strain': phi = phi_0 + eps_1, "
       "the strain-based update of Liu et al. (2024, CMAME 429:117165) eq. (40), where eps_1 "
       "is the maximum principal strain of `strain_property`.");
   params.addParam<MaterialPropertyName>(
       "strain_property", "mechanical_strain",
       "Rank-two kinematic strain tensor whose maximum principal value drives the strain-based "
       "porosity update. Only used when porosity_update_model = strain. Defaults to the "
       "'mechanical_strain' property declared by ComputeSmallStrain.");
   ```
   - Update the existing `addClassDescription(...)` text to mention both laws.
   - `phase_field` stays `addRequiredCoupledVar` (unchanged). It is unused by the strain
     branch but kept required to avoid an interface change; all current callers already pass it.

3. **Constructor changes**
   - Initialize the new members in the initializer list:
     ```cpp
     _porosity_update_model(getParam<MooseEnum>("porosity_update_model") == "strain"
                                ? PorosityUpdateModel::STRAIN
                                : PorosityUpdateModel::DAMAGE),
     _mechanical_strain(nullptr)
     ```
     (Place these after the existing members; `_porosity_damaged` stays last to preserve the
     existing declaration that depends on `isNodal()`.)
   - In the constructor body, keep the existing `lower > upper` check, then add:
     ```cpp
     if (_porosity_update_model == PorosityUpdateModel::STRAIN)
     {
       if (isNodal())
         mooseError("porosity_update_model = strain is only available at quadrature points; "
                    "the qp strain tensor is unavailable for the nodal porosity property in ",
                    name(), ". Use porosity_update_model = damage for nodal porosity.");
       _mechanical_strain = &getMaterialProperty<RankTwoTensor>("strain_property");
     }
     ```
     `getMaterialProperty` must be called only in this branch so the damage model does not
     create a dependency on `mechanical_strain` (keeping the default path usable in inputs that
     have no strain material).

4. **`computeQpProperties()` rewrite** (preserve the final clamp shared by both branches):
   ```cpp
   Real phi_raw;
   if (_porosity_update_model == PorosityUpdateModel::STRAIN)
   {
     std::vector<Real> eigvals(LIBMESH_DIM);
     (*_mechanical_strain)[_qp].symmetricEigenvalues(eigvals);
     const Real eps1 = eigvals.back();            // ascending order -> last is max principal strain
     phi_raw = _initial_porosity + eps1;          // Liu et al. (2024) eq. (40)
   }
   else
   {
     const Real dmg = std::clamp(_damage[_qp], 0.0, 1.0);
     const Real g = std::pow(1.0 - dmg, 2);       // valid for AT1/AT2
     phi_raw = _initial_porosity + (1.0 - _initial_porosity) * (1.0 - g);
   }
   _porosity_damaged[_qp] = std::clamp(phi_raw, _porosity_lower_bound, _porosity_upper_bound);
   ```

### Interfaces
- New input parameters on `ElkPorousFlowDamagedPorosity`:
  `porosity_update_model = damage|strain` (default `damage`),
  `strain_property = <MaterialPropertyName>` (default `mechanical_strain`).
- No change to declared output property names or their types.

### Edge Cases to Handle
- `porosity_update_model = strain` with `isNodal()` → `mooseError` (see 3).
- `porosity_update_model = strain` but `mechanical_strain` (or the named `strain_property`)
  not declared by any material → MOOSE raises the standard "material property not found"
  error via `getMaterialProperty`; acceptance test must confirm a clear failure (no silent
  zero strain). No extra code needed beyond calling `getMaterialProperty`.
- Pure compression / 2D plane strain: `eps_1 = max(..., 0) >= 0` so `phi_raw >= phi_0`;
  combined with `lower_bound` the porosity never drops below `min(phi_0, ...)`. Verify in tests.
- `phi_raw` exceeding `upper_bound` (large opening) → clamped to `upper_bound`. Verify.
- `damage` branch must remain byte-identical to current output (regression test).

### Acceptance Criteria
- [ ] `farms-opt` compiles with the modified material (or `make` succeeds where the toolchain
      links; see Risk note about the macOS Tahoe link issue — compilation of the object file is
      the binding criterion if linking is blocked locally).
- [ ] An input with no `porosity_update_model` parameter reproduces current damaged-porosity
      values exactly (Phase 3 `damage_model` test).
- [ ] With `porosity_update_model = strain` and a prescribed uniform `eps_xx = a > 0`
      (`eps_yy = eps_zz = 0`), the qp porosity equals `clamp(phi_0 + a, lower, upper)`
      (Phase 3 `strain_model` test).
- [ ] `porosity_update_model = strain` with `isNodal()` produces the documented `mooseError`.

### Dependencies
- Depends on: nothing.
- Required by: Phase 2 (AD mirror), Phase 3 (tests), Phase 4 (docs/usage).

---

## Phase 2: Mirror the option in the AD twin `ElkADPorousFlowDamagedPorosity` — DEFERRED

> **Deferred for this pass** (per the locked scope decision). Do not implement now. Documented
> here so the future AD mirror matches Phase 1 exactly. Skip straight to Phase 3 after Phase 1.

### Goal
Keep the AD material feature-identical so AD-based assemblies can also select the strain law.

### Files to Modify
- `include/materials/porousflowmatprops/ElkADPorousFlowDamagedPorosity.h`
- `src/materials/porousflowmatprops/ElkADPorousFlowDamagedPorosity.C`

### Detailed Requirements
1. Header: add `#include "RankTwoTensor.h"`; add the same `PorosityUpdateModel` enum and
   `const PorosityUpdateModel _porosity_update_model;` plus
   `const ADMaterialProperty<RankTwoTensor> * _mechanical_strain;`.
2. `validParams()`: add the identical `porosity_update_model` and `strain_property` params
   (copy text from Phase 1).
3. Constructor: parse the enum identically; in the `STRAIN` branch guard `isNodal()` and bind
   `_mechanical_strain = &getADMaterialProperty<RankTwoTensor>("strain_property");`.
4. `computeQpProperties()`: branch as in Phase 1 but with AD types:
   ```cpp
   ADReal phi_raw;
   if (_porosity_update_model == PorosityUpdateModel::STRAIN)
   {
     std::vector<ADReal> eigvals(LIBMESH_DIM);
     (*_mechanical_strain)[_qp].symmetricEigenvalues(eigvals);
     const ADReal eps1 = eigvals.back();
     phi_raw = _initial_porosity + eps1;
   }
   else
   {
     ADReal dmg = std::min(std::max(_damage[_qp], ADReal(0.0)), ADReal(1.0));
     ADReal g = pow(1.0 - dmg, 2);
     phi_raw = _initial_porosity + (1.0 - _initial_porosity) * (1.0 - g);
   }
   _porosity_damaged[_qp] = std::min(std::max(phi_raw, ADReal(_porosity_lower_bound)),
                                     ADReal(_porosity_upper_bound));
   ```

### Edge Cases to Handle
- The AD strain branch requires an **AD** `mechanical_strain` provider (e.g.
  `ADComputeSmallStrain`). Document this in the class comment and param help: pairing the AD
  material's strain branch with a non-AD `ComputeSmallStrain` will fail the
  `getADMaterialProperty` lookup. The `damage` branch has no such requirement (unchanged).

### Acceptance Criteria
- [ ] AD material compiles.
- [ ] AD `damage` branch reproduces the existing
      `pulsepower/unit_tests/phase1_materials/test_ad_damaged_porosity.i` values exactly.
- [ ] AD `strain` branch, fed an AD strain with known `eps_1`, returns `clamp(phi_0 + eps_1)`.

### Dependencies
- Depends on: Phase 1 (for the agreed enum/param names and semantics).
- Required by: nothing. **This phase is secondary** — production inputs use the non-AD
  material; implement only after Phase 1 + Phase 3 pass. It may be deferred if no AD consumer
  needs the strain law.

---

## Phase 3: Unit tests (formal regression suite)

### Goal
Lock in the unchanged damage behavior and verify eq. (40) analytically, including the clamps.

### Files to Create
- `test/tests/materials/damaged_porosity/damage_model.i` — non-AD regression of the damage law.
- `test/tests/materials/damaged_porosity/strain_model.i` — analytic strain-law verification.
- `test/tests/materials/damaged_porosity/tests` — MOOSE test spec.
- `test/tests/materials/damaged_porosity/gold/` — gold files (generated by running the inputs
  once the code is in and values are confirmed against the hand calculations below).
- `test/tests/materials/damaged_porosity/EXPECTED_VALUES.md` — written derivation of every
  gold number (mirror the style of
  `test/tests/materials/permeability_normal_strain/EXPECTED_VALUES.md`).
- (Optional) `doc/content/source/materials/ElkPorousFlowDamagedPorosity.md` — short design page
  so the `design =` field in `tests` resolves. If a quicker path is preferred, point `design`
  at an existing materials doc page; do not leave `design` pointing at a non-existent file.

### Detailed Requirements

1. **`damage_model.i`** — mirror `pulsepower/unit_tests/phase1_materials/test_ad_damaged_porosity.i`
   but with the **non-AD** `ElkPorousFlowDamagedPorosity` and a `MaterialRealAux` (not
   `ADMaterialRealAux`). Prescribe `d = x` via `FunctionIC`/`FunctionAux`. Omit
   `porosity_update_model` (exercise the default). Use `phi_0 = 0.008`,
   `porosity_lower_bound = 0.008`, `porosity_upper_bound = 0.999`. Sample porosity at
   `x = 0.1667, 0.5, 0.8333` with `PointValue` postprocessors. Expected (hand-computed):
   - `x≈0.1667`: `g=(1-0.1667)^2=0.6944`, `phi=0.008+0.992*0.3056=0.3112`.
   - `x=0.5`:   `g=0.25`, `phi=0.008+0.992*0.75=0.752`.
   - `x≈0.8333`: `g=(1-0.8333)^2=0.02779`, `phi=0.008+0.992*0.97221=0.97444`.
   Output `csv = true`; test with `CSVDiff`.

2. **`strain_model.i`** — impose a known strain without solving elasticity, following the
   pattern in `test/tests/materials/permeability_normal_strain/isotropic_axis.i`:
   - `[GlobalParams] displacements = 'disp_x disp_y'`.
   - `disp_x`, `disp_y` as `AuxVariables` (FIRST LAGRANGE) set by `FunctionAux`.
   - `[Materials]`: `ComputeSmallStrain` (declares `mechanical_strain`), and the
     `ElkPorousFlowDamagedPorosity` block with `porosity_update_model = strain`,
     `phase_field = d` (d can be a constant aux = 0), `initial_porosity = 0.008`,
     `porosity_lower_bound = 0.008`, `porosity_upper_bound = 0.999`. No elasticity/stress
     material is required — only the strain tensor is consumed.
   - A trivial `[Variables] u` with a `Diffusion` kernel + Dirichlet BCs to make the `Steady`
     solve non-empty (exactly as the permeability tests do).
   - Mesh: `GeneratedMeshGenerator dim=2 nx=3 ny=1 xmin=0 xmax=1`. Impose a **piecewise-linear**
     `disp_x` whose slope (= `eps_xx`) is constant within each element by placing breakpoints on
     element boundaries `x = 0, 1/3, 2/3, 1`. Using `PiecewiseLinear` with
     `x = '0 0.33333333 0.66666667 1'` and cumulative `y` values chosen so element slopes are:
     - element 0 (center x=1/6): slope `1e-3`  → `eps_xx = 1e-3`,
     - element 1 (center x=1/2): slope `2.0`   → `eps_xx = 2.0` (drives the upper clamp),
     - element 2 (center x=5/6): slope `-0.5`  → `eps_xx = -0.5` (compression).
     `disp_y = 0` so `eps_yy = eps_zz = 0`.
   - Extract porosity into a `MONOMIAL CONSTANT` aux via `MaterialRealAux`
     (`property = PorousFlow_porosity_qp_damaged`) and sample with `PointValue` at element
     centers `(1/6, 0.5, 0)`, `(1/2, 0.5, 0)`, `(5/6, 0.5, 0)`.
   - Expected (eq. 40, then clamp), eigenvalues of `diag(eps_xx, 0, 0)`:
     - elem 0: `eps_1 = max(1e-3, 0, 0) = 1e-3` → `phi = 0.008 + 1e-3 = 0.009`.
     - elem 1: `eps_1 = 2.0` → `phi_raw = 2.008` → clamp upper → `0.999`.
     - elem 2: `eps_1 = max(-0.5, 0, 0) = 0` → `phi = 0.008` (compression → initial porosity).
   - Output `csv = true`; test with `CSVDiff`, `rel_err = 1e-8`.

3. **Lower-bound clamp** — add one more test entry reusing `strain_model.i` with
   `cli_args = 'Materials/<porosity_block>/porosity_lower_bound=0.05'` and a distinct
   `Outputs/file_base`, asserting elem 0 porosity clamps up to `0.05`
   (`phi_raw = 0.009 < 0.05`). Provide its own gold CSV. (If `cli_args` overriding nested
   block names is awkward, instead add a second tiny input `strain_lower_clamp.i`.)

4. **Nodal-guard error test** — add a `RunException` test entry: an input requesting the strain
   law on a nodal porosity instance, expecting the `isNodal()` `mooseError`. If constructing a
   nodal instance in a standalone input is impractical, document this as covered by the missing-
   strain `getMaterialProperty` failure instead and drop this entry. Do not fabricate a passing
   test.

5. **`tests` spec** — one `[Tests]` block with `design = '<doc page>.md'`,
   `issues = '#strain_based_porosity'` (or the tracking id the team uses), and `requirement`
   strings for each entry, mirroring `test/tests/materials/permeability_normal_strain/tests`.

### Edge Cases to Handle
- Element-center `PointValue` on a `MONOMIAL CONSTANT` field must sample the element actually
  containing the point; keep `nx=3` so the three sample points land in distinct elements and
  strains are element-uniform (linear disp).
- Ensure `eps_zz` truly is 0 in 2D (`ComputeSmallStrain` plane strain) — the compression case
  relies on it for `eps_1 = 0`.

### Acceptance Criteria
- [ ] `./run_tests --re damaged_porosity` passes all entries.
- [ ] `damage_model` gold matches the hand calculations in step 1.
- [ ] `strain_model` gold matches `0.009`, `0.999`, `0.008` in step 2.
- [ ] Lower-bound entry yields `0.05` for element 0.
- [ ] `EXPECTED_VALUES.md` derives every gold number.

### Dependencies
- Depends on: Phase 1 (non-AD) for `strain_model`/`damage_model`; Phase 2 only if an AD test is
  added (not required).
- Required by: Phase 4.

---

## Phase 4: Documentation and comparison usage

### Goal
Make the bounded-vs-strain comparison runnable by changing a single parameter, and document it.

### Files to Modify / Create
- `include/.../ElkPorousFlowDamagedPorosity.h` (and AD header): expand the class doxygen block
  to describe both laws and cite eq. (40).
- `doc/content/source/materials/ElkPorousFlowDamagedPorosity.md` (if created in Phase 3):
  document parameters, the two laws, the qp-only constraint, the AD-strain requirement, and a
  worked example of switching between bounded-maximum and strain-based porosity.
- One existing production elasticity input (e.g.
  `pulsepower/cmame_revision/2d_hydromech/permeability_formula/elasticity_E1d25.i`): add the new
  parameter to the `[porosity_damaged]` block **as a commented option** so the default behavior
  is untouched, e.g.:
  ```
  [porosity_damaged]
    type = ElkPorousFlowDamagedPorosity
    phase_field = d
    initial_porosity = ${porosity}
    porosity_lower_bound = 0.008
    porosity_upper_bound = 0.999
    # porosity_update_model = strain   # eq. (40): phi = phi0 + max principal strain
    #                                   # (default 'damage' = bounded-maximum porosity)
  []
  ```

### Detailed Requirements
1. Do **not** alter any default-run behavior of production inputs (keep the new line commented).
2. In the docs, state explicitly that the comparison study is: run the case once with
   `porosity_update_model = damage` (bounded-maximum baseline) and once with
   `porosity_update_model = strain` (eq. 40), holding everything else fixed, then compare the
   `porosity_aux` field and any porosity-dependent outputs (Biot modulus, fluid kinetic/elastic
   energy).
3. Note that `mechanical_strain` is the kinematic (total) strain in the absence of thermal or
   plastic eigenstrains; if a thermal eigenstrain is later added, switch `strain_property` to
   `total_strain` to keep eq. (40) representing the kinematic crack opening.

### Acceptance Criteria
- [ ] Production inputs run unchanged with the commented option present.
- [ ] Docs explain how to reproduce the bounded-vs-strain comparison.

### Dependencies
- Depends on: Phase 1 (and Phase 3 if the doc page backs the `design` field).
- Required by: nothing.

---

## Testing Strategy
- **Phase 1/2:** compile the app; rely on Phase 3 for behavioral verification.
- **Phase 3 regression (`damage_model`):** guarantees the default path is byte-stable, so no
  existing production run changes.
- **Phase 3 analytic (`strain_model`):** verifies eq. (40) and both clamps against
  hand-computed values for prescribed, element-uniform strains — no dependence on an elasticity
  solve, so the numbers are exact.
- **Error-path test:** the `isNodal()` guard (or, failing that, the missing-`mechanical_strain`
  failure) is exercised so misuse fails loudly.
- Run `./run_tests --re damaged_porosity` and, for the regression guarantee, the existing
  `permeability_normal_strain` and `phase1_materials` suites to confirm no collateral change.

## Risk Assessment
- **Nodal porosity + strain law.** PorousFlow can request nodal porosity; the qp `mechanical_strain`
  is unavailable there. Mitigation: hard `mooseError` in the constructor (Phase 1 step 3).
  Detection: the nodal-guard test. Production inputs only use the qp instance, so this does not
  affect current runs.
- **Jacobian coupling (non-AD).** The non-AD material returns a plain `Real` with no strain
  derivative, so in a monolithic HM Newton solve the strain-based porosity lags within an
  iteration (consistent with how the damage-based porosity is already treated explicitly via the
  staggered phase-field subapp). Risk: slightly slower nonlinear convergence in strongly coupled
  steps. Detection: monitor nonlinear iteration counts. Mitigation if needed: use the AD twin
  (Phase 2) with an AD strain provider for full coupling.
- **Strain choice (`mechanical_strain` vs `total`/`elastic`).** For the current HM model with no
  thermal eigenstrain these coincide; the configurable `strain_property` param lets the user
  retarget if that changes. Documented in Phase 4.
- **Reversibility.** Eq. (40) is instantaneous: if a crack closes (strain relaxes) porosity
  decreases, unlike the monotonic damage law (damage is irreversible). This is physically
  intended and is the effect the study aims to observe; if a monotone (max-over-history) variant
  is later wanted, it is a clean extension (add a stateful `Old` property and take
  `max(phi_old, phi_new)`) — out of scope here, noted for awareness.
- **Eigenvalue ordering.** The plan assumes `symmetricEigenvalues` returns ascending order
  (verified: it delegates to LAPACK `dsyev`). The `strain_model` upper/compression cases would
  fail loudly if this assumption were wrong, so the test guards it.
- **Toolchain link issue (local).** Per project memory, `farms-opt` may fail to *link* in the
  conda MOOSE env on macOS 26.x though objects still compile. Treat successful compilation of the
  modified translation units as the local binding criterion; run the `run_tests` suite wherever
  the executable links (CI / Linux box) for the behavioral acceptance criteria.
```
