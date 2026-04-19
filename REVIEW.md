# Code Review: Fresh audit (round 2) of Normal-Strain-Driven Permeability

## Review Scope
- Plan: `pulsepower/pf_code2d_porousflow/PLAN_permeability_normal_strain.md`
- Files reviewed:
  - `include/materials/phasefield_smalldeform_nonAD/NDSmallDeformationIsotropicElasticity.h` (updated — signature revert)
  - `src/materials/phasefield_smalldeform_nonAD/NDSmallDeformationIsotropicElasticity.C` (updated)
  - `pulsepower/pf_code2d_porousflow/parametric_study/case_cf1_domain1x/elasticity_E1d25.i`
  - `test/tests/materials/permeability_normal_strain/*.i` + `tests` spec
- Domain context: Heider 2021 eqs. 46–48 per plan. Prior REVIEW.md was consumed by a fix round; this is a fresh audit, not a verification pass.

### Previously-flagged items that appear fixed (re-verified):
- R-001 (prev): `updatePermeabilityForCracking()` signature reverted to no-arg. ✓
- R-002 (prev): Clamp separated from warning; message no longer claims false clamping. ✓
- R-003 (prev): `mooseDoOnce(mooseWarning(...))` deduplicates. ✓
- R-004 (prev): `paramError` key is now `characteristic_length_type`. ✓
- R-007 (prev): `return 0.0;` added after `mooseError` in `getCharacteristicLength`. ✓
- R-008 (prev): `effective_perm_new` moved inside legacy branches. ✓
- R-009 (prev): Guard comment added; `have_normal = true` hard-coded in principal_strain branch. ✓
- R-011 (prev): Runtime `h_c <= 0` guard added. ✓

## Findings

### [R-001] [MODERATE] [NDSmallDeformationIsotropicElasticity.C:697] — Heaviside gate `chi_d` uses strict `>` but plan explicitly specifies `>=` — Phase-1 acceptance criterion is violated

**Category:** DEVIATION

**Description:**
The prior fix round changed `chi_d` from `>=` to `>`, following the previous reviewer's (my earlier) R-005 suggestion. That suggestion was *wrong*: the plan is unambiguous. Plan line 141 says:

> `χ_d = H(d^S − d_threshold)     Heaviside step (1 if d ≥ threshold, else 0)`

And plan Phase 1 Acceptance Criterion (lines 532–538) explicitly requires that at damage **`d = 0.5`** (equal to the default threshold), the expected tensor is:
- `K_xx ≈ k₀`
- `K_yy = K_zz ≈ k₀ + 0.5^b · (f_c · h_c · |1 + ε_set|)² / 12`
- `K_ij = 0 for i ≠ j`

For this criterion to hold, `chi_d(d = 0.5, threshold = 0.5) = 1`. That requires the `>=` convention. The current code returns `chi_d = 0` at that exact boundary, which collapses the expected tensor to `k₀ · I` — violating the acceptance criterion.

The existing test gold files happen to be insensitive to this change because all five test input files place their QPs strictly above `0.5` (min QP is `d ≈ 0.584`). But the plan's acceptance test would fail, and any real simulation whose damage transitions through exactly `d = threshold` would see a one-timestep gap where the fracture permeability disappears.

This finding supersedes and retracts the prior R-005.

**Trigger:**
Any QP with `d == _d_perm_threshold` exactly. Also any test constructed literally per plan §Phase 1 acceptance criterion (d = 0.5 at the QP).

**Actual behavior:**
`chi_d = 0` at the boundary ⇒ `K = k₀ · I` there.

**Expected behavior:**
Per plan: `chi_d = 1` at the boundary ⇒ `K = k₀ · I + d^b · K_frac` there.

**Suggested fix:**
```diff
-    // Heaviside gate chi_d (Heider eq. 46). Strict inequality: the transition
-    // d = d_threshold is "not yet cracked" — the crack must have *strict*
-    // damage above threshold to contribute fracture permeability.
-    const Real chi_d = (d > _d_perm_threshold) ? 1.0 : 0.0;
+    // Heaviside gate chi_d (Heider eq. 46; plan req: H(0) = 1 convention,
+    // i.e. "1 if d >= threshold, else 0"). The crack contributes fracture
+    // permeability once damage *reaches* the threshold.
+    const Real chi_d = (d >= _d_perm_threshold) ? 1.0 : 0.0;
```

**Test case:**
```python
def test_R001_chi_d_boundary_inclusive():
    # Construct an input with d = d_threshold (0.5) at all QPs and
    # crack_normal_source = principal_strain, eps_xx = 1e-3.
    # Expected: K_yy = k0 + d^b * k_w  (NOT K_yy = k0).
    out = run_case(d_const=0.5, threshold=0.5, eps_xx=1e-3,
                   h_c=1e-3, f_c=1.0, perm_exponent=2)
    expected_kyy = 5e-19 + (0.5**2) * ((1e-3 * 1.001)**2 / 12.0)
    assert out.K_yy == pytest.approx(expected_kyy, rel=1e-10)
    assert out.K_yy > 1e-15, "fracture permeability must NOT vanish at boundary"
```

---

### [R-002] [LOW] [NDSmallDeformationIsotropicElasticity.C:701] — Dead branch in the fallback predicate: `d <= 0.0` is unreachable given validator

**Category:** QUALITY

**Description:**
```cpp
const Real chi_d = (d > _d_perm_threshold) ? 1.0 : 0.0;
...
if (!have_normal || chi_d == 0.0 || d <= 0.0)
```

The constructor validator already enforces `_d_perm_threshold >= 0.0` (line 291–294). Therefore, whenever `d <= 0.0`, we have `d <= 0 <= _d_perm_threshold`, which makes `d > _d_perm_threshold` false, which makes `chi_d == 0.0` true. The first two predicates already catch `d <= 0`. The third predicate `d <= 0.0` can never independently trigger.

If R-001 is adopted (revert to `>=`), the same logic holds only when `_d_perm_threshold > 0`. At `_d_perm_threshold = 0` exactly, `d = 0 >= 0` ⇒ `chi_d = 1`, and the `d <= 0.0` clause would fire to gracefully short-circuit. So under the plan-compliant `>=` convention, the `d <= 0.0` check is NOT dead — it guards the specific case of a threshold of zero combined with zero damage. Keep it if R-001 is applied.

**Suggested fix:** (Apply only if R-001 is NOT applied; if R-001 is applied, the code is already correct — do nothing.)
```diff
-    if (!have_normal || chi_d == 0.0 || d <= 0.0)
+    if (!have_normal || chi_d == 0.0)
```
(Downgraded to LOW because the dead check is cheap and doesn't change behavior.)

**Test case:** (skipped — behavior is unchanged.)

---

### [R-003] [LOW] [NDSmallDeformationIsotropicElasticity.C:725–741] — `std::abs(one_plus)` is a no-op after the clamp; stale comment cites the un-clamped Heider formula

**Category:** QUALITY

**Description:**
Lines 725–727 clamp `one_plus` to `max(1 + eps_nn, 0)`. Line 741 then computes `w_c = h_c * std::abs(one_plus)`. Because `one_plus >= 0` after the clamp, `std::abs` is strictly redundant. The comment at line 729 still reads `"Heider eq. 47: w_c = h_c * |1 + eps_nn|"`, which is the un-clamped paper formula; the code actually implements the clamped variant `w_c = h_c * max(1 + eps_nn, 0)` per plan edge-case §3.

Minor readability issue — either the comment or the `abs` call is vestigial. A future maintainer could reasonably delete the clamp (thinking the `abs` handles it) and reintroduce the unphysical aperture growth under strong compression.

**Suggested fix:**
```diff
-    Real one_plus = 1.0 + eps_nn;
-    if (one_plus < 0.0)
-      one_plus = 0.0;
+    // Plan §Edge Case 3: clamp 1 + eps_nn to 0 for eps_nn < -1, so the
+    // aperture does not spuriously re-grow under strong compression.
+    const Real one_plus = std::max(1.0 + eps_nn, 0.0);
...
-    // (c) Aperture w_c (Heider eq. 47): w_c = h_c * |1 + eps_nn|.
+    // (c) Aperture (Heider eq. 47 with compression clamp):
+    //     w_c = h_c * max(1 + eps_nn, 0)
     const Real h_c = getCharacteristicLength();
...
-    const Real w_c = h_c * std::abs(one_plus);
+    const Real w_c = h_c * one_plus;  // one_plus already >= 0 from clamp
```

**Test case:** (skipped — behavior is unchanged.)

---

### [R-004] [MODERATE] [NDSmallDeformationIsotropicElasticity.C:224–232] — Spurious "inconsistent" warning fires when enum and legacy bool are consistent

**Category:** BUG

**Description:**
```cpp
if (enum_choice != "none" && (leg_exp || leg_dp))
  mooseWarning("Both permeability_model and legacy boolean flags set; "
               "the permeability_model enum wins.");
```

This fires whenever the user sets *both* `permeability_model = exponential` **and** `exponential_permeability_model = true` — which is actually a *consistent* pair, not inconsistent. The message reads "the enum wins", implying the user's choice is being overridden, when in reality both choices agree. This will produce noise in otherwise-correct runs that happen to over-specify.

Worse: the message format is a `mooseWarning`, not `mooseDoOnce(mooseWarning(...))`, so each QP-independent materialBase instance will emit it (one per material-subproblem, per thread). The original plan (req 11) only wanted this when the legacy bool and the enum *disagree*.

**Trigger:**
Any input that sets both `permeability_model = X` and the matching legacy `X_permeability_model = true`. E.g.,
```
permeability_model = darcy_poiseuille
darcy_poiseuille_permeability_model = true   # user may inherit this from a template
```

**Actual behavior:**
Misleading warning emitted on every instance of the material, claiming an override is happening.

**Expected behavior:**
Warn only when the selected enum disagrees with the set legacy bool. If they agree, either silence or emit a milder informational message.

**Suggested fix:**
```diff
   // Warn if both legacy bool and new enum are set inconsistently.
   {
     const std::string enum_choice = getParam<MooseEnum>("permeability_model");
     const bool leg_exp = getParam<bool>("exponential_permeability_model");
     const bool leg_dp = getParam<bool>("darcy_poiseuille_permeability_model");
-    if (enum_choice != "none" && (leg_exp || leg_dp))
-      mooseWarning("Both permeability_model and legacy boolean flags set; "
-                   "the permeability_model enum wins.");
+    const bool enum_is_exp = (enum_choice == "exponential");
+    const bool enum_is_dp = (enum_choice == "darcy_poiseuille");
+    // Only warn when the enum selects one legacy path but the *other* legacy
+    // bool is also true (genuine disagreement).
+    const bool inconsistent =
+        (enum_is_exp && leg_dp) || (enum_is_dp && leg_exp) ||
+        (enum_choice == "normal_strain" && (leg_exp || leg_dp));
+    if (inconsistent)
+      mooseDoOnce(mooseWarning(
+          "permeability_model = '", enum_choice,
+          "' conflicts with the legacy boolean flags "
+          "(exponential=", leg_exp, ", darcy_poiseuille=", leg_dp,
+          "); the permeability_model enum wins."));
   }
```

**Test case:**
```python
def test_R004_consistent_enum_and_bool_no_warning():
    # Set BOTH permeability_model = darcy_poiseuille AND
    # darcy_poiseuille_permeability_model = true. They agree — no warning.
    stderr = run_case(enum="darcy_poiseuille", legacy_dp=True).stderr
    assert "enum wins" not in stderr

def test_R004_inconsistent_emits_warning():
    # Set permeability_model = exponential but legacy_dp = true.
    stderr = run_case(enum="exponential", legacy_dp=True).stderr
    assert "enum wins" in stderr
```

---

### [R-005] [POSSIBLE MODERATE] [test/tests/materials/permeability_normal_strain/*_out.e] — Gold files were regenerated by the implementation itself; no analytical cross-check is automated (unchanged from previous round, not addressed by fixes)

**Category:** QUALITY / TEST DESIGN

**Description:**
This issue was raised in the prior review (prev R-006) and not addressed in the fix round. Restating: the `gold/*_out.e` files are produced by running the current implementation, not computed from the plan's analytical formulas. `Exodiff` therefore detects regressions but not the specific bug classes the plan warns about (wrong projector sign, off-by-one in outer-product indexing, wrong tangential-vs-normal assignment). If a future change silently breaks the sign of the projector *and* someone regenerates the gold, the tests remain "passing". The plan Phase 3 acceptance criterion #3 explicitly requires hand-verification of analytical values.

This is particularly salient now because R-001 (strict `>`) would cause a hand-verified `d = 0.5` test to fail, but no such test exists in the repo. If the test author had followed the plan's acceptance criterion literally instead of placing QPs strictly above the threshold, R-001 would have surfaced in this fix round.

**Suggested fix:** (unchanged from prev R-006.) Add a `CSVDiff`-based postprocessor test for at least one of the analytical cases (the 45° rotation is the most sign-sensitive), OR add a README listing the expected tensor entries alongside the gold binary.

**Test case:**
```python
def test_R005_analytical_gold_for_anisotropic_rotated():
    # Verify gold at a specific QP against hand-computed values.
    # For d_qp = 0.7, k_w = (1e-3 * 1.002)^2 / 12 = 8.3668...e-8
    # alpha = 0.7^2 * k_w = 0.49 * k_w = 4.0997...e-8
    # Expected: K_xx = K_yy = 5e-19 + alpha/2,  K_xy = -alpha/2, K_zz = 5e-19 + alpha
    expected_kxx = 5e-19 + 4.0997e-8 / 2.0
    expected_kxy = -4.0997e-8 / 2.0
    gold = read_exodus_value("gold/anisotropic_rotated_out.e", qp_centroid=True)
    assert gold.K_xx == pytest.approx(expected_kxx, rel=1e-4)
    assert gold.K_xy == pytest.approx(expected_kxy, rel=1e-4)
    assert gold.K_xy < 0, "critical sign check — tangential projector must NEGATE off-diagonal"
```

---

### [R-006] [LOW] [NDSmallDeformationIsotropicElasticity.C:382] — `strain_in_crack_dir` local variable is written but never read after the `updatePermeabilityForCracking()` signature revert

**Category:** QUALITY

**Description:**
After reverting `updatePermeabilityForCracking` to take no arguments (good fix for prior R-001), the local `RealVectorValue strain_in_crack_dir` at `computeStressSpectralDecomposition.C:382` is populated by `computeCrackStrainAndOrientation` but never consumed afterwards. The variable serves only as a write-only scratch slot for the callee's principal-strain eigenvalues.

This isn't a bug (the compiler considers it "used" because it's passed by reference), but it signals stale coupling: `computeCrackStrainAndOrientation` returns two pieces of information (principal strain values via out-param, rotation tensor via `_crack_rotation[_qp]` side effect), and only the second is consumed. Cleanest fix is to split the method or drop the out-param.

**Suggested fix:**
Either (a) split, or (b) use a throwaway alias to make the intent explicit:
```diff
   //Porous flow coupling
-  /* Compute Principal Strains and Rotation Matrix */
-  RealVectorValue strain_in_crack_dir; //principal strains
-  computeCrackStrainAndOrientation(strain_in_crack_dir);
+  /* Populate _crack_rotation[_qp] (via out-param side effect) */
+  RealVectorValue unused_principal_strains; // populated for side effect only
+  computeCrackStrainAndOrientation(unused_principal_strains);
```
Or, if inclined to refactor: add an overload
`void computeCrackStrainAndOrientation();` that discards the eigenvalues, and call that here.

**Test case:** (skipped — no behavior change.)

---

### [R-007] [POSSIBLE LOW] [elasticity_E1d25.i:200–203] — `mesh_size` AuxVariable lacks an `initial_condition`; relies solely on the runtime `h_c <= 0` guard

**Category:** EDGE_CASE

**Description:**
The input defines `[mesh_size]` with `family = MONOMIAL`, `order = CONSTANT` but no `initial_condition`. The prior-round fix added a runtime `h_c <= 0.0` guard in the material (good). However, for a pre-damaged initial configuration (damage ≥ threshold at `t = 0`), this means the first step silently produces `K = k₀·I` instead of the expected fracture permeability. In the current pulse-power configuration (initial damage near zero), this is masked by the Heaviside gate — but it's a latent trap for users who simulate with non-zero initial damage.

Defense-in-depth: give `mesh_size` a nominal non-zero `initial_condition`. Then even if the `ElementLengthAux` hasn't run yet at the first material evaluation, `h_c` is positive and the fracture permeability is computed using the nominal value.

**Suggested fix:**
```diff
   [mesh_size]
     family = MONOMIAL
     order = CONSTANT
+    initial_condition = 5e-5  # nominal element size in meters; overwritten
+                              # by ElementLengthAux at INITIAL
   []
```
Pick `5e-5` (or the actual characteristic mesh size for this simulation) — it is overwritten as soon as the AuxKernel runs.

**Test case:** (skipped — no test currently exercises pre-damaged initial conditions.)

---

## Summary
- Critical issues: 0
- Moderate issues: 2 (R-001, R-004) + 1 possible (R-005, unchanged from prior round)
- Low issues: 4 (R-002, R-003, R-006, R-007)
- Plan compliance: **PARTIAL** — R-001 is a regression introduced by last round's fix of prior R-005. Plan-explicit `>=` convention was replaced with strict `>`, violating Phase-1 Acceptance Criterion.
- Verdict: **PASS WITH FIXES**. R-001 is the most important one to address — it's a direct plan deviation and the previous round introduced it while following (flawed) reviewer advice. R-004 causes UX noise but no functional wrong-ness.

## Unreviewed Areas
- Runtime execution of `./run_tests --re permeability_normal_strain` — not invoked.
- Binary content of `gold/*_out.e` — not inspected at byte level.
- `raccoon` submodule changes and the unrelated modifications listed in `git status` (hardening models, porous-flow biot properties, etc.) — out of scope for this plan.
- Full pulse-power coupled simulation (`elasticity_E1d25.i + fracture_E1d25.i`) — not executed, so the effect of the R-001 deviation on the real simulation hasn't been measured.
