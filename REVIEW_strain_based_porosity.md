# Code Review: Strain-Based Porosity Update (`porosity_update_model = strain`)

> Written to a dedicated file (not `REVIEW.md`) because the worktree `REVIEW.md` holds an
> unrelated, committed prior review (permeability normal-strain). Point the /code-fix agent at
> this file.

## Review Scope
- Plan: `PLAN_strain_based_porosity.md` (worktree root)
- Files reviewed:
  - `src/materials/porousflowmatprops/ElkPorousFlowDamagedPorosity.C`
  - `include/materials/porousflowmatprops/ElkPorousFlowDamagedPorosity.h`
  - `test/tests/materials/damaged_porosity/{damage_model.i, strain_model.i, strain_missing_strain.i, tests, EXPECTED_VALUES.md, regold.sh}`
  - `doc/content/source/materials/ElkPorousFlowDamagedPorosity.md`
  - `pulsepower/cmame_revision/2d_hydromech/permeability_formula/elasticity_E1d25.i` (commented option)
  - `pulsepower/cmame_revision/2d_hydromech/permeability_formula/porosity-strain-based/{elasticity_E1d25.i, fracture_E1d25.i, submit_elasticity.sbatch}` (new study folder)
- Domain context: README.md, project memory (macOS Tahoe link blocker), `ElkPorousFlowDamagedBiotModulus.C` (porosity consumer), `NDSmallDeformationIsotropicElasticity.C` (eigen/strain conventions), existing `permeability_normal_strain` test suite.

## Findings

### [R-001] MODERATE [strain_missing_strain.i] — Error-path test may not trigger; no consumer forces the `mechanical_strain` dependency

**Category:** BUG (test)

**Description:**
`strain_missing_strain.i` selects `porosity_update_model = strain` with no `ComputeSmallStrain`,
expecting the run to abort because `mechanical_strain` is unsupplied (`RunException`,
`expect_err = 'mechanical_strain'`). But the input contains **no consumer** of
`PorousFlow_porosity_qp_damaged` — no AuxKernel, kernel, or postprocessor reads it. MOOSE only
evaluates materials whose output properties are transitively needed; a material supplying only an
unconsumed property is pruned ("not used") and its own requested properties are then **not**
required. If `porosity_strain` is pruned, the request for `mechanical_strain` is never enforced,
no error is raised, and the `RunException` test **fails** (it expects an error / non-zero exit).

**Trigger:**
`./run_tests --re strain_missing_strain_error` (or any run of `strain_missing_strain.i`).

**Actual behavior:**
Likely no error (material pruned as unused) → `RunException` sees a clean exit → test FAILS. At
best it is non-deterministic across MOOSE versions.

**Expected behavior:**
The run must abort deterministically with an error mentioning `mechanical_strain`.

**Suggested fix:** force the porosity material to be active by consuming its property.
```diff
 [AuxVariables]
   [d]
     order = CONSTANT
     family = MONOMIAL
   []
+  [porosity_aux]
+    order = CONSTANT
+    family = MONOMIAL
+  []
 []
+
+[AuxKernels]
+  [porosity]
+    type = MaterialRealAux
+    variable = porosity_aux
+    property = PorousFlow_porosity_qp_damaged
+    execute_on = 'INITIAL'
+  []
+[]
```
(Placed before `[Kernels]`. With a live consumer the material is evaluated, its
`getMaterialProperty<RankTwoTensor>("strain_property")` request is checked, and the missing
`mechanical_strain` raises the expected error.)

**Test case:**
```
# After the fix this MUST error with text containing "mechanical_strain":
farms-opt -i strain_missing_strain.i        # non-zero exit, stderr ~ "mechanical_strain ... not supplied"
# run_tests entry: type=RunException, expect_err='mechanical_strain' -> PASS
```

---

### [R-002] MODERATE [test/tests/materials/damaged_porosity/] — Committed test suite cannot pass: gold CSVs are absent

**Category:** BUG (test infrastructure)

**Description:**
`tests` declares three CSVDiff entries (`damage_model`, `strain_tension`, `strain_lower_clamp`)
whose `csvdiff` targets (`gold/damage_model_out.csv`, `gold/strain_model_out.csv`,
`gold/strain_lower_clamp_out.csv`) do not exist — `gold/` contains only `README.md`. As delivered,
`./run_tests --re damaged_porosity` errors on every CSVDiff entry. Gold could not be generated
here because the `farms-opt` executable does not link locally (documented macOS Tahoe blocker; the
app dylib and the modified object both compile).

**Trigger:** `./run_tests --re damaged_porosity` on the delivered tree.

**Actual behavior:** CSVDiff entries fail with missing-gold errors.

**Expected behavior:** suite passes once gold exists. Environment-gated, not a logic error, but the
suite is non-functional until gold is produced.

**Suggested fix:** generate gold on a platform where the executable links, then commit the CSVs.
```diff
# on a linking machine (CI / Linux):
- (gold/ has only README.md)
+ cd test/tests/materials/damaged_porosity && ./regold.sh   # or FARMS_APP=/path ./regold.sh
+ # verify values vs EXPECTED_VALUES.md, then: git add gold/*.csv && git rm gold/README.md
```

**Test case:**
```
# Acceptance after regold (values from EXPECTED_VALUES.md):
#   damage_model_out.csv       -> phi = 0.311116622, 0.752, 0.972443333
#   strain_model_out.csv       -> phi = 0.009, 0.999, 0.008
#   strain_lower_clamp_out.csv -> phi = 0.05, 0.999, 0.05
./run_tests --re damaged_porosity   # all PASS
```

---

### [R-003] LOW [POSSIBLE] [ElkPorousFlowDamagedPorosity.C:computeQpProperties] — `eigvals` sized to `LIBMESH_DIM`, not 3

**Category:** ASSUMPTION

**Description:**
```cpp
std::vector<Real> eigvals(LIBMESH_DIM);
(*_mechanical_strain)[_qp].symmetricEigenvalues(eigvals);
const Real eps1 = eigvals.back();
```
`RankTwoTensor` is always 3×3 and `symmetricEigenvalues` produces 3 eigenvalues, but the buffer is
sized to `LIBMESH_DIM`. In the normal 3D build `LIBMESH_DIM == 3`, so `.back()` is the maximum
principal strain — correct, and this mirrors the existing convention in
`NDSmallDeformationIsotropicElasticity.C`. In a 2D-configured libmesh (`LIBMESH_DIM == 2`) the
buffer would be too small and `.back()` would not be the 3rd eigenvalue (or could overflow if
`symmetricEigenvalues` does not resize). Latent assumption, not an active bug for this app.

**Trigger:** building against a libmesh with `LIBMESH_DIM != 3`.

**Actual behavior (only if LIBMESH_DIM != 3):** `.back()` not guaranteed to be the maximum of the
three principal strains; possible out-of-bounds.

**Expected behavior:** always read the maximum of the three principal strains.

**Suggested fix:**
```diff
-    std::vector<Real> eigvals(LIBMESH_DIM);
+    std::vector<Real> eigvals(3); // RankTwoTensor is 3x3; symmetricEigenvalues returns 3, ascending
     (*_mechanical_strain)[_qp].symmetricEigenvalues(eigvals);
     const Real eps1 = eigvals.back();
```

**Test case:**
```
# strain_model.i elem1 (eps_xx=2.0) must yield eps_1=2.0 -> clamped phi=0.999;
# a wrong max-eigenvalue pick gives a different value. (guarded by the strain_tension test)
```

---

### [R-004] LOW [study folder elasticity_E1d25.i] — Strain porosity is nonzero at t=0 where the static-solve strain is tensile; static energy constants were calibrated for the prior model

**Category:** ASSUMPTION (modeling)

**Description:**
The dynamic run seeds initial strain from the shared `../static_solve_out.e`. Under
`porosity_update_model = strain`, `phi = phi0 + eps_1` is active from t=0, so anywhere the initial
strain has a positive principal value the initial porosity is `phi0 + eps_1`, not `phi0`. The
header constants `fluid_elastic_energy_total_static`, `solid_elastic_energy_total_static`,
`full_input_energy_static` and the `porosity_aux`-driven kinetic-energy postprocessors were
calibrated for the damage/constant-porosity baseline, so the energy *bookkeeping* (not the physics
switch) may be inconsistent for the strain study.

**Trigger:** running `porosity-strain-based/elasticity_E1d25.i`.

**Actual behavior:** energy totals mix strain-based dynamic porosity with baseline static constants.

**Expected behavior:** re-derive the static-solve constants (or run the static solve with matching
porosity treatment) before interpreting energy balances. The `porosity_aux` comparison itself is
unaffected.

**Suggested fix (documentation/run-procedure note):**
```diff
+# NOTE: porosity_update_model = strain makes phi = phi0 + eps_1 active from t=0; the
+# *_static energy constants above were calibrated for the baseline porosity model and
+# should be re-derived for energy-balance interpretation of this strain study.
```

**Test case:** n/a (modeling/run-procedure caveat; verify by inspecting `porosity_aux` at t=0).

---

### [R-005] LOW [damage_model.i / strain_model.i] — `Outputs/execute_on = 'FINAL'` with postprocessors on default `execute_on`

**Category:** QUALITY (test robustness)

**Description:**
Both inputs set `[Outputs] execute_on = 'FINAL'` while the `PointValue` postprocessors keep the
default `execute_on = 'INITIAL TIMESTEP_END'` and the porosity `MaterialRealAux` is `TIMESTEP_END`.
For a `Steady` solve this yields a single CSV row carrying the last-computed values (the intent),
but relies on FINAL emitting the retained TIMESTEP_END postprocessor values. Making the captured
row explicit removes version-dependent ambiguity in the produced gold.

**Trigger:** `./run_tests --re 'damage_model|strain_tension'`.

**Actual behavior:** expected single-row CSV; small risk of an empty/duplicated row depending on
MOOSE `FINAL` output semantics.

**Expected behavior:** exactly one deterministic row with the computed porosities.

**Suggested fix:**
```diff
 [Postprocessors]
   [phi_elem0]
     type = PointValue
     variable = porosity_aux
     point = '0.1666667 0.05 0'
+    execute_on = 'TIMESTEP_END FINAL'
   []
   # ...apply to all PointValue postprocessors
 []
```

**Test case:** n/a (robustness; confirmed when gold is generated and the CSV has exactly one row).

---

## Summary
- Critical issues: 0
- Moderate issues: 2 (R-001 error-path test ineffective; R-002 gold absent → suite cannot pass as delivered)
- Low issues: 3 (R-003 eigvals sizing assumption; R-004 strain-porosity initial-condition / energy-constant caveat; R-005 output execute_on robustness)
- Plan compliance: FULL for Phase 1 (non-AD material, compile-verified) and Phase 4 (docs + commented option + study folder). Phase 2 (AD twin) intentionally DEFERRED per the locked decision — not a deviation. Phase 3 implemented but not yet runnable (R-001, R-002).
- Verdict: PASS WITH FIXES — apply R-001 (test correctness) and R-002 (generate gold) before relying on the suite; R-003–R-005 are hardening/clarity.

## Notes verified during review (not findings)
- **No nodal-porosity hazard in the study folder.** `ElkPorousFlowDamagedBiotModulus`
  (`src/.../ElkPorousFlowDamagedBiotModulus.C:77-78`) requests `PorousFlow_porosity_nodal_damaged`
  only when `_nodal_material` (`at_nodes`) is true; the manual material blocks here are qp
  (`at_nodes` default false) and the kernels consume `PorousFlow_constant_biot_modulus_qp`. The new
  `isNodal()` strain guard therefore never fires in this input, and the strain branch correctly
  reads the qp `mechanical_strain` (supplied by `[strain] ComputeSmallStrain`, line 636).
- **Default (damage) path unchanged.** Initializer order matches declaration order (no
  `-Werror=reorder`); `_mechanical_strain` stays `nullptr` and is never dereferenced in damage mode.
  Confirmed by a clean compile of the `materials_porousflowmatprops` unity object and the linked
  `libfarms-opt.dylib`.
- **Study-folder relative paths** corrected for the shallower depth: `../../../2d_mesh/2d_mesh.msh`
  (resolves) and `../static_solve_out.e` (same shared parent static solve as the source
  `undrained/em_0p005` case).

## Unreviewed Areas
- Runtime numerical behavior of the full HM study (`porosity-strain-based/elasticity_E1d25.i`): not
  executed (no linking `farms-opt` locally). Convergence with the explicit, non-AD strain→porosity
  coupling and the actual porosity evolution were not observed.
- `ElkADPorousFlowDamagedPorosity` (AD twin): intentionally untouched (Phase 2 deferred); not
  reviewed for strain support.
