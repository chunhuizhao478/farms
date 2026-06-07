# Code Review: regularize_crack_normal feature + strained-based run (2026-06-05)

## Review Scope
- Plan: conversational `/code-implement` prompt + user clarifications in-thread
  (Case 0 "smooth everywhere", d=0 vs d=1 separated by the chi_d gate, and the
  `strained-based` production run using normal-regularization with eps=500).
  No formal PLAN.md for this change.
- Files reviewed:
  - `src/materials/phasefield_smalldeform_nonAD/NDSmallDeformationIsotropicElasticity.C`
  - `include/materials/phasefield_smalldeform_nonAD/NDSmallDeformationIsotropicElasticity.h`
  - `test/tests/materials/permeability_normal_strain/regularized_normal_core.i` (new)
  - `test/tests/materials/permeability_normal_strain/tests` (entry added)
  - `test/tests/materials/permeability_normal_strain/EXPECTED_VALUES.md` (Test 8 added)
  - `pulsepower/cmame_revision/2d_hydromech/permeability_formula/strained-based/`
    (`elasticity_E1d25.i`, `fracture_E1d25.i`, `submit_elasticity.sbatch`)
- Domain context: Heider (2021) eqs. 46-48 normal-strain permeability, prior
  `EXPECTED_VALUES.md`, the in-repo macOS-Tahoe toolchain note (executable link
  is broken; tests could not be run / gold could not be generated this session).

NOTE: this file supersedes the previous REVIEW.md (residual-aperture review). If
the residual_aperture R-001 (non-portable gold) was never addressed, re-check it
separately - it is out of scope for this change.

## Findings

### [R-001] [MODERATE] [tests / regularized_normal_core] — Exodiff test references a gold file that does not exist

**Category:** BUG (test infrastructure)

**Description:**
`test/tests/materials/permeability_normal_strain/tests` adds an `Exodiff` test
`[regularized_normal_core]` with `exodiff = 'regularized_normal_core_out.e'`, but
there is no `gold/regularized_normal_core_out.e`. Confirmed: `ls gold/ | grep
regular` returns nothing. An Exodiff test with no gold file errors at run time
("gold file does not exist"). The implementer documented that the gold could not
be generated because the executable will not link on this macOS (Tahoe 26.5)
toolchain, but the test entry was still committed pointing at a nonexistent gold,
so `./run_tests` will report this test as errored for everyone until the gold is
produced.

**Trigger:**
`./run_tests -i permeability_normal_strain` (once the app builds) -> the
`regularized_normal_core` test errors on the missing gold.

**Actual behavior:**
Test errors (no gold), counted as a failure in the suite.

**Expected behavior:**
Either the gold file exists and the test passes, or the test is explicitly
deferred so it does not error.

**Suggested fix:**
Preferred: once the app links, generate the gold and verify against Test 8 in
`EXPECTED_VALUES.md` (`K_xx=K_yy=K_zz≈4.083333e-8`, off-diag 0):
```bash
cd test/tests/materials/permeability_normal_strain
farms-opt -i regularized_normal_core.i && cp regularized_normal_core_out.e gold/
```
If it must be committed before the toolchain is fixed, gate it so it does not
error in the meantime - add to the test block in `tests`:
```diff
  [regularized_normal_core]
    type = 'Exodiff'
    input = 'regularized_normal_core.i'
    exodiff = 'regularized_normal_core_out.e'
+   # TODO(remove once gold generated): app executable cannot link on macOS
+   # Tahoe 26.5 toolchain, so the gold could not be produced this session.
+   skip = 'gold pending: regularized_normal_core_out.e not yet generated'
    requirement = '...'
```

**Test case:**
```
# After generating the gold, this must PASS:
./run_tests -i permeability_normal_strain --re regularized_normal_core
# Assert: status OK, eff_perm_00 == eff_perm_11 == eff_perm_22 ≈ 4.083333e-8,
#         eff_perm_01 ≈ 0  (isotropic despite permeability_anisotropic = true).
```

---

### [R-002] [MODERATE] [tests] — No test for the d < threshold (undamaged) separation, which the user explicitly required

**Category:** EDGE_CASE / missing test

**Description:**
The whole point of the user's clarification ("separate d=0 with d=1") is that the
regularized normal must NOT enhance permeability in undamaged material: at
`d < damage_threshold_for_permeability` with `grad(d)=0` the regularized path
sets `n_d=0, have_normal=true`, and correctness then depends *entirely* on the
downstream `chi_d == 0` gate routing the point to `k0*I`. The only new test
(`regularized_normal_core.i`) uses `d = 0.7` (above threshold -> isotropic
fracture perm). There is no test exercising the undamaged branch with
`regularize_crack_normal = true`. A future regression that reordered/removed the
`chi_d` gate, or that let `have_normal=true` bypass it, would silently turn the
undamaged bulk into fracture perm and no test would catch it.

**Trigger:**
Any change to the `if (!have_normal || chi_d == 0.0 || d <= 0.0)` gate ordering
or condition.

**Actual behavior:**
Undamaged-side separation is unverified.

**Expected behavior:**
A test asserting that `regularize_crack_normal = true` + uniform `d < threshold`
yields `K = k0*I` (matrix perm, no enhancement).

**Suggested fix:**
Add `test/tests/materials/permeability_normal_strain/regularized_normal_undamaged.i`
(clone `regularized_normal_core.i`, change only `d_aux` function `'0.7'` -> `'0.3'`),
and a `tests` entry. Expected: `eff_perm_00=eff_perm_11=eff_perm_22 = 5e-19`
(= k0), off-diag 0. Add the new case to the EXPECTED_VALUES regeneration loop.
```diff
  [d_aux]
    type = FunctionAux
    variable = d
-   function = '0.7'
+   function = '0.3'   # below damage_threshold_for_permeability (0.5)
    execute_on = 'INITIAL TIMESTEP_END'
  []
```

**Test case:**
```
# regularized_normal_undamaged.i (uniform d=0.3, regularize_crack_normal=true)
# Assert: eff_perm_00 == eff_perm_11 == eff_perm_22 == 5e-19 (k0), eff_perm_01==0.
# i.e. the regularized normal does NOT create fracture perm below threshold.
```

---

### [R-003] [MODERATE] [POSSIBLE] [validParams / constructor] — Default crack_normal_regularization = 1e-8 makes the feature a silent no-op for realistic gradients

**Category:** ASSUMPTION / usability bug

**Description:**
`crack_normal_regularization` defaults to `1e-8` (1/length). For phase-field
fracture, `|grad(d)| ~ 1/l`; in the production run `l = 2e-4 m` => `|grad(d)| ~
5e3 /m` across the damage band. The shrink factor `s = |grad d|/(|grad d| + eps)`
is then `5e3/(5e3 + 1e-8) ≈ 1 - 2e-12` - i.e. `n_d` is the unit normal to ~12
digits everywhere except exactly at `grad(d)=0`. At finite-element quadrature
points the gradient is essentially never `<= 1e-8`, so a user who enables
`regularize_crack_normal = true` and leaves `crack_normal_regularization` at the
default gets the legacy anisotropic behavior with no visible isotropic core - the
feature appears enabled but does nothing. The docstring even says "Choose it
small relative to |grad(d)| ~ 1/l," which, taken with the `1e-8` default, points
users toward an ineffective value. (The production input correctly overrides it
to 500; the trap is for anyone relying on the default.)

**Trigger:**
`regularize_crack_normal = true` with `crack_normal_regularization` left at the
default and any realistic `|grad(d)| >> 1e-8`.

**Actual behavior:**
Regularization silently inactive; isotropic-core behavior never appears.

**Expected behavior:**
Enabling the feature should either do something visible, or fail loudly so the
user knows they must pick `eps` on the gradient scale.

**Suggested fix:**
Require `crack_normal_regularization` to be set explicitly when the feature is on,
in the `normal_strain` constructor block (next to the existing regularize check):
```diff
     if (_regularize_crack_normal &&
         _normal_source != CrackNormalSource::damage_gradient)
       paramError("regularize_crack_normal",
                  "regularize_crack_normal = true only applies to "
                  "crack_normal_source = damage_gradient. ...");
+    if (_regularize_crack_normal &&
+        !isParamSetByUser("crack_normal_regularization"))
+      paramError("crack_normal_regularization",
+                 "regularize_crack_normal = true requires "
+                 "crack_normal_regularization to be set explicitly. The default "
+                 "(1e-8 /m) is far below the phase-field gradient scale "
+                 "|grad(d)| ~ 1/l, so the regularization would be a no-op. Set "
+                 "eps on the order of the near-core |grad(d)| you want to "
+                 "isotropize (e.g. a fraction of 1/l).");
```
(Alternatively keep the default but emit a one-time `mooseWarning` when
`regularize_crack_normal=true` and the param was not set by the user.)

**Test case:**
```
# Input with: permeability_model=normal_strain, porous_flow_coupling=true,
#   crack_normal_source=damage_gradient, regularize_crack_normal=true,
#   and NO crack_normal_regularization line.
# Assert: construction fails with a paramError naming crack_normal_regularization
#   (or, for the warning variant, a mooseWarning is emitted).
```

---

### [R-004] [LOW] [updatePermeabilityForCracking] — Stale "unit crack normal" comment and redundant d<=0 gate for the regularized path

**Category:** QUALITY

**Description:**
Two minor issues introduced/exposed by the change, neither alters results:
1. Line ~736 comment "(a) Determine unit crack normal n_d." is no longer accurate
   for the regularized branch: `n_d = grad(d)/(|grad(d)|+eps)` has magnitude
   `|grad d|/(|grad d|+eps) < 1`, i.e. it is deliberately NOT a unit vector. A
   reader checking "is n_d unit?" will be misled.
2. In the fallback `if (!have_normal || chi_d == 0.0 || d <= 0.0)`, the `d <= 0.0`
   term is dead for the normal-strain path because `chi_d == 0.0` already fires
   for every `d < _d_perm_threshold` (default 0.5), which includes all `d <= 0`.
   Harmless but invites confusion about whether `d` can be negative here.

**Trigger:** N/A (readability only).

**Actual behavior:** Comment overstates that n_d is unit; redundant condition.

**Expected behavior:** Comment reflects that n_d may be sub-unit under
regularization.

**Suggested fix:**
```diff
-    // (a) Determine unit crack normal n_d.
+    // (a) Determine crack normal n_d. NOTE: with regularize_crack_normal = true
+    //     this is grad(d)/(|grad(d)|+eps), whose magnitude is < 1 (it tends to 0
+    //     at the crack core); it is a true unit vector only on the legacy paths.
```
(The redundant `d <= 0.0` can be left as defensive, or dropped - low priority.)

**Test case:** N/A (comment-only).

---

### [R-005] [LOW] [POSSIBLE] [updatePermeabilityForCracking] — At the core the aperture loses all normal-strain dependence

**Category:** ASSUMPTION (design consequence, user-accepted)

**Description:**
Because `eps_nn = n_d·eps·n_d` uses the regularized `n_d`, as `n_d -> 0` at the
core `eps_nn -> 0`, so `w_c = h_c·|1+eps_nn| -> h_c` regardless of the actual
strain state. The "normal-strain permeability model" therefore degenerates at the
crack core to a strain-independent isotropic perm `K = k0*I + d^b·(h_c^2/12)·I` -
a highly-strained open core and a barely-open core get the same aperture `h_c`.
This is the direct consequence of the Case-0 design the user explicitly approved,
so it is recorded, not asserted as a bug. Flagging so it is a conscious choice and
not later mistaken for a regression. If physical aperture at the core is desired,
the model would need `w_c` driven by something other than `n_d·eps·n_d` when
`|n_d|` is shrunk (e.g. use the unshrunk unit normal for `eps_nn` but the shrunk
normal only for the projector) - a larger design change, not recommended now.

**Trigger:** Fully-damaged core (`grad(d)->0`) under non-zero strain.

**Actual/Expected:** Accepted behavior; no change requested.

**Suggested fix:** None (documentation/awareness only). Optionally note this in
the validParams docstring for `regularize_crack_normal`.

**Test case:** N/A.

---

## Summary
- Critical issues: 0
- Moderate issues: 3 (R-001 missing gold, R-002 missing undamaged-side test,
  R-003 ineffective default eps)
- Low issues: 2 (R-004 comment/redundancy, R-005 design note)
- Plan compliance: PARTIAL — source logic matches the approved Case-0 design and
  the d=0/d=1 separation is correct; test coverage is incomplete (no gold, no
  undamaged-side test) and the default-eps ergonomics undercut the feature.
- Verdict: PASS WITH FIXES — the production source path is correct; fix R-001/
  R-002/R-003 before relying on the test suite or shipping the feature defaults.

## Unreviewed Areas
- Runtime/numerical verification: the executable cannot be linked on this macOS
  Tahoe 26.5 toolchain, so no test was actually run and no gold was generated.
  All correctness claims here (and the EXPECTED_VALUES Test 8 numbers) are from
  static reading + hand calc, not execution.
- `strained-based/` run: relative paths and the regularization block were checked
  statically and resolve correctly; the actual solver behavior / eps=500 efficacy
  was not run (and depends on this run's damage profile - see R-003 rationale).
- Prior REVIEW.md content (residual_aperture) was superseded by this file; its
  open items, if any, were not re-verified here.
