# Code Review: porosity_bounded worktree (2026-06-05)

> Written as `POROSITY_BOUNDS_REVIEW.md` (not `REVIEW.md`) to avoid clobbering
> the repo's existing `REVIEW.md`, which belongs to the unrelated "Degraded Bulk
> Modulus" effort — same rationale as `POROSITY_BOUNDS_PLAN.md` vs `PLAN.md`.

## Review Scope
- Plan: `POROSITY_BOUNDS_PLAN.md` (worktree root)
- Files reviewed:
  - `test/tests/materials/damaged_porosity/POROSITY_BOUNDS_REFERENCE.md` (new)
  - `POROSITY_BOUNDS_PLAN.md` (new)
  - `src/materials/porousflowmatprops/ElkPorousFlowDamagedPorosity.C` (description strings)
  - `src/materials/porousflowmatprops/ElkADPorousFlowDamagedPorosity.C` (description strings)
  - `pulsepower/unit_tests/phase1_materials/test_ad_damaged_porosity.i` (comment)
  - `pulsepower/cmame_revision/2d_hydromech/permeability_formula/porosity_bounded/{elasticity_E1d25.i, fracture_E1d25.i, submit_elasticity.sbatch}` (new case)
- Domain context: companion `permeability_normal_strain/EXPECTED_VALUES.md`; sibling
  cases `permeability_formula/undrained/em_0p005` (source) and `em_0p010` (depth/path control).

## Verified non-issues (highest-risk items, checked explicitly)
- **Relative-path conversions in the new case are correct.** Source `em_0p005`
  is 3 levels below `2d_hydromech`; new `porosity_bounded` is 1 level below, so
  each path drops two `../`. Confirmed against the canonical sibling `em_0p010`
  (`mesh = ../../static_solve_out.e` from depth-3) — the new
  `mesh = ../permeability_formula/static_solve_out.e` (from depth-1) resolves to
  the *same* `2d_hydromech/permeability_formula/static_solve_out.e`, produced by
  `permeability_formula/static_solve.i`. `../../2d_mesh/2d_mesh.msh` resolves to
  the existing `cmame_revision/2d_mesh/2d_mesh.msh`. The only "missing" target,
  `static_solve_out.e`, is a generated artifact absent in every case until the
  static solve runs — not a defect.
- **No stray `0.999` left in the new case** (grep clean); the only
  `porosity_upper_bound` is the intended `0.065`.
- **Bound ordering valid:** lower `0.008` ≤ upper `0.065`, and
  `initial_porosity = ${porosity} = 0.008 ∈ [0.008, 0.065]`, so the material
  ctor's `lower > upper` guard does not fire and no clamp pathology arises.
  Kinetic-energy postprocessor `…/(porosity_aux²)` stays finite (φ ≥ 0.008).
- **C++ edits are string-literal-only.** Adjacent-literal concatenation produces
  correctly-spaced help text; defaults (`0.0`, `0.999`), range checks, and the
  `std::clamp`/`std::min/std::max` logic are untouched.
- **Test comment edit is cosmetic** — `test_ad_damaged_porosity.i` samples
  d ≤ 0.833 (φ ≤ 0.972), never reaching the clamp; results unaffected.

## Findings

### [R-001] MODERATE `test/tests/materials/damaged_porosity/POROSITY_BOUNDS_REFERENCE.md` — reference doc contradicts the implemented 0.065 case

**Category:** DEVIATION

**Description:**
The reference document is the project's authoritative justification for the
porosity cap, and it argues — repeatedly and explicitly — to **keep `0.999`**:
"Bottom line up front. `0.999` is … the numerical regularization of the model's
`phi -> 1` void limit"; §4 "My recommendation: keep `0.999`"; and it argues
*against* a low physical cap: "lowering the cap is … physically inconsistent with
treating the crack as an open conduit for flow." But the new production input
`porosity_bounded/elasticity_E1d25.i` implements the opposite decision,
`porosity_upper_bound = 0.065` (the measured-max physical cap from §3). The
repo's justification doc and its implemented case now tell contradictory
stories; a reader cannot tell which is authoritative.

**Trigger:** Anyone reading `POROSITY_BOUNDS_REFERENCE.md` alongside the
`porosity_bounded` case.

**Actual behavior:** Doc recommends `0.999` and argues against low caps; case uses `0.065`.

**Expected behavior:** The doc must record the `porosity_bounded` decision
(0.065 = Staněk & Géraud measured maximum, used as a physically-bounded
skeleton-cap sensitivity case) and reconcile its recommendation, e.g. "default
0.999 = void-limit regularization; the `porosity_bounded` case deliberately
caps at 0.065 to model the damaged zone as a still-granular solid — accepting
the storage-vs-conduit inconsistency noted above as the point of that test."

**Suggested fix:** Append a subsection to §4 (or a new §4a) and amend the
bottom-line:
```diff
- My recommendation: keep `0.999` and justify it as the void-limit regularization
+ Default: `0.999` (void-limit regularization). The `porosity_bounded` sensitivity
+ case (`pulsepower/cmame_revision/2d_hydromech/permeability_formula/porosity_bounded/`) instead caps at
+ `0.065` — the maximum *measured* fractured/altered granite porosity (§3, Staněk &
+ Géraud 2019) — to model the damaged zone as a still-granular solid skeleton rather
+ than an open void. This deliberately trades the conduit interpretation (§4 caveat)
+ for a physically-bounded porosity; compare the two runs.
```

**Test case (verification):**
```bash
# Must mention the implemented case + value, and must not leave a bare "keep 0.999":
grep -q "porosity_bounded" test/tests/materials/damaged_porosity/POROSITY_BOUNDS_REFERENCE.md
grep -q "0.065"            test/tests/materials/damaged_porosity/POROSITY_BOUNDS_REFERENCE.md
! grep -qi "recommendation: keep .0.999" test/tests/materials/damaged_porosity/POROSITY_BOUNDS_REFERENCE.md
```

---

### [R-002] MODERATE `POROSITY_BOUNDS_PLAN.md` — completion record omits the delivered `porosity_bounded` case

**Category:** DEVIATION

**Description:**
The plan/completion record marks every phase "✅ complete" but documents only the
reference doc + C++ description edits. It does not mention the most substantial
delivered artifact: the `porosity_bounded` case folder, the `0.065` cap decision,
the relative-path conversions, or the kept cluster base path. The record
therefore overstates "all complete" while omitting the work the user actually
requested last, and a future reader/`/code-fix` consumer gets an inaccurate
picture of branch contents.

**Trigger:** Reading `POROSITY_BOUNDS_PLAN.md` to understand what the branch contains.

**Actual behavior:** No reference to `porosity_bounded`, `0.065`, or path conversions.

**Expected behavior:** A phase/section recording the new case and its decisions.

**Suggested fix:** Append:
```diff
+ ## Phase 4 — porosity_bounded sensitivity case  ✅
+ - Copied `permeability_formula/undrained/em_0p005/{elasticity_E1d25.i,
+   fracture_E1d25.i, submit_elasticity.sbatch}` → `2d_hydromech/permeability_formula/porosity_bounded/`.
+ - Set `porosity_upper_bound = 0.065` (measured-max physical cap; see reference doc §3).
+ - Corrected relative paths for the shallower folder depth (drop two `../`):
+   mesh `../../../../2d_mesh/2d_mesh.msh` → `../../2d_mesh/2d_mesh.msh` (both inputs);
+   solution `../../static_solve_out.e` → `../permeability_formula/static_solve_out.e`.
+ - Submit script input path repointed to `porosity_bounded`; cluster base
+   `/scratch2/.../farms_cdms_04192026/` kept (merge to cluster later).
+ - OPEN: submit job name still `cmame_2dhm_undr_0p005` (see review R-003).
```

**Test case (verification):**
```bash
grep -q "porosity_bounded" POROSITY_BOUNDS_PLAN.md
grep -q "0.065"            POROSITY_BOUNDS_PLAN.md
```

---

### [R-003] LOW `pulsepower/cmame_revision/2d_hydromech/permeability_formula/porosity_bounded/submit_elasticity.sbatch` — job/log names still identify the source case

**Category:** QUALITY

**Description:**
The copied submit script keeps `-J cmame_2dhm_undr_0p005` and output/error
patterns `cmame_2dhm_undr_0p005.o%j`/`.e%j`. The job now runs the
`porosity_bounded` (0.065) case, so the queue name and log files misidentify it
and (if launched from a shared CWD) interleave with the original `em_0p005`
case's logs. The `%j` job-id suffix prevents hard filename collision, but the
naming is misleading. (User is aware; left per "only the path".)

**Trigger:** Submitting this script; inspecting `squeue`/log files.

**Actual behavior:** Job/log names say `undr_0p005`.

**Expected behavior:** Names reflect the new case.

**Suggested fix:**
```diff
-#SBATCH -J cmame_2dhm_undr_0p005            # Job name
-#SBATCH -o cmame_2dhm_undr_0p005.o%j        # Name of stdout output file
-#SBATCH -e cmame_2dhm_undr_0p005.e%j        # Name of stderr error file
+#SBATCH -J cmame_2dhm_poro_0p065            # Job name
+#SBATCH -o cmame_2dhm_poro_0p065.o%j        # Name of stdout output file
+#SBATCH -e cmame_2dhm_poro_0p065.e%j        # Name of stderr error file
```

---

### [R-004] LOW [POSSIBLE] `porosity_bounded/elasticity_E1d25.i` — 0.065 cap is physically inconsistent with the aperture-based cubic-law permeability

**Category:** ASSUMPTION

**Description:**
The damaged permeability is a crack-aperture cubic law (`k = w²/12`, open
conduit) while `porosity_upper_bound = 0.065` keeps ~93.5 % solid fraction in the
fully-damaged zone for the poroelastic storage / Biot modulus. This is the exact
storage-vs-conduit inconsistency the reference doc warns about. It is a
deliberate user modeling choice (sensitivity test), **not a code defect** — flagged
per review rule 7 so it is not silently lost: the run's damaged-zone storativity
will not correspond to an open crack. Reasoning that it is *possible* rather than
confirmed-bug: the cap may be intentionally exploring a granular-skeleton regime.

**Trigger:** Running the `porosity_bounded` case and interpreting damaged-zone
fluid storage/pressure.

**Suggested fix:** No code change required. Ensure R-001's doc update records the
intent so the inconsistency is understood as the purpose of the case, not an
oversight.

---

### [R-005] LOW `ElkPorousFlowDamagedPorosity.C` / `ElkADPorousFlowDamagedPorosity.C` — help text hardcodes "0.999"

**Category:** QUALITY

**Description:**
The `porosity_upper_bound` description reads "…; 0.999 regularizes the
fully-damaged phi->1 void limit…". It is accurate for the *default*, but now that
a production input overrides it to `0.065`, a user who set a different value sees
help text asserting `0.999`. Minor; reframing to "the default 0.999" removes the
ambiguity.

**Trigger:** `--show-input` / `--dump` on an input that overrides the bound.

**Suggested fix (both files, identical):**
```diff
-      "Upper clamp applied to the updated porosity; 0.999 regularizes the "
+      "Upper clamp applied to the updated porosity; the default 0.999 regularizes the "
       "fully-damaged phi->1 void limit (see "
```

---

## Summary
- Critical issues: 0
- Moderate issues: 2 (R-001 doc/impl contradiction, R-002 stale plan record)
- Low issues: 3 (R-003 job name, R-004 modeling-consistency note, R-005 help text)
- Plan compliance: PARTIAL — the original doc+code tasks are FULL, but the plan
  record (R-002) does not cover the later-delivered `porosity_bounded` case, and
  the reference doc (R-001) contradicts it.
- Verdict: PASS WITH FIXES — no functional/path defects in the new case; fix the
  two MODERATE documentation-consistency findings before merge so the branch's
  justification matches its implementation.

## Unreviewed Areas
- Runtime correctness of the `porosity_bounded` case (actual MOOSE solve) — not
  executed; no build tree in this worktree and `static_solve_out.e` is generated
  on the cluster. Static path/parameter checks only.
- The cluster-side deployment (whether `porosity_bounded/` is synced to
  `/scratch2/.../farms_cdms_04192026/`) — out of scope; user will merge later.
