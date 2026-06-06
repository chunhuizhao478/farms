# Plan & completion record: justify `porosity_upper_bound` in `ElkPorousFlowDamagedPorosity`

> Worktree `porosity-bounded` (branch `porosity-bounded`, from `cdms`).
> This file is kept separate from the repo root `PLAN.md`, which belongs to the
> unrelated "Degraded Bulk Modulus" effort and must not be modified here.

## Goal
Justify the numerically-convenient `porosity_upper_bound = 0.999` (and the
`porosity_lower_bound = 0.008`) of `ElkPorousFlowDamagedPorosity` against the
literature, and document the justification in the codebase.

## Deliverables (all complete)
- `test/tests/materials/damaged_porosity/POROSITY_BOUNDS_REFERENCE.md` — cited
  justification (new).
- Doc-pointer comments in the two material sources' `validParams()`
  (`ElkPorousFlowDamagedPorosity.C`, `ElkADPorousFlowDamagedPorosity.C`) —
  description strings only; defaults/logic unchanged.
- Tightened a misleading comment in
  `pulsepower/unit_tests/phase1_materials/test_ad_damaged_porosity.i`
  (`phi_raw = 1.0 -> clamped to 0.999`, was loosely `= 1.0 (clamped...)`).

## Phase 1 — Codebase audit  ✅
- Both materials default `lower=0.0`, `upper=0.999`
  (`ElkPorousFlowDamagedPorosity.C:17-20`, `ElkADPorousFlowDamagedPorosity.C:26-29`).
- Production: 89 `.i` set `upper=0.999`; 83 set `lower=0.008` (4 use `${porosity}`=0.008).
- Outliers (non-production, recorded in the reference doc §5): two archived
  sensitivity inputs use `upper=0.9`/`0.5`; three kernel test fixtures use
  `lower=0.1`.
- `lower_bound (0.008) == initial_porosity ${porosity} (0.008)` confirmed in production.
- Downstream `phi -> 1` singularities identified: damaged Biot modulus
  `denom = phi/K_f + (alpha-phi)/K_grain` (negative as `phi -> 1`, `alpha < 1`;
  `ElkPorousFlowDamagedBiotModulus.C:115-122`) and `(1-phi)` / `1/phi^2`
  solid-fraction terms (`ElkADThreeFieldHistoryEnergyEnhanced.C`).

## Phase 2 — Literature review  ✅
- Intact granite: Schild et al. (2001) GJI 146, 111–125 — unaltered 0.8–1.2 vol%,
  in-situ 0.5%, lab 1.0–1.17%. Westerly granite ≈1.4% (corroborating).
  → `phi_0 = 0.008` (0.8%) at low end of intact range.
- Damaged granite: Staněk & Géraud (2019) Solid Earth 10, 251–274 — fractured/
  altered max measured porosity 6.5% (≪ 0.999).
- Phase-field framework: Heider (2021) Eng. Fract. Mech. 253, 107881 — the repo's
  own permeability-model reference; fully-cracked zone = fluid-filled void →
  `phi -> 1` limit, regularized by the 0.999 cap (retains 1e-3 solid fraction).
- All recorded with DOIs and VERIFIED/UNVERIFIED tags in the reference doc §6.

## Phase 3 — Document in the codebase  ✅
- Reference doc finalized (6 sections: model & limits, lower bound, "what 0.999
  is not", upper bound regularization, codebase audit, references).
- Identical description edits applied to both AD and non-AD `validParams()`;
  `git diff` shows only string literals changed (no default, range check, or
  clamp logic touched).
- Explicit honesty: doc states 0.999 is NOT a measured granite porosity but the
  regularized model void limit.

## Phase 4 — `porosity_bounded` sensitivity case  ✅
- Copied `permeability_formula/undrained/em_0p005/{elasticity_E1d25.i,
  fracture_E1d25.i, submit_elasticity.sbatch}` →
  `2d_hydromech/permeability_formula/porosity_bounded/`.
- Set `porosity_upper_bound = 0.065` (measured-max physical cap; reference doc §3/§4a).
- Relative paths set for the final folder depth (1 level below `permeability_formula`,
  matching sibling `oil_water_mixture/`):
  mesh `../../../2d_mesh/2d_mesh.msh` (both inputs);
  solution `../static_solve_out.e` (shared parent static-solve output, as in `em_0p005`).
- Submit script `-i` path points to `permeability_formula/porosity_bounded/`; cluster
  base `/scratch2/.../farms_cdms_04192026/` kept (will merge to cluster later); job/log
  names renamed to `cmame_2dhm_poro_0p065`.

## Verification notes / caveats
- The material/doc edits (Phases 1–3) are description-string + Markdown only → no
  behavioral change; the existing unit test's sampled points (`d ≤ 0.833`,
  `phi ≤ 0.972`) never reach the clamp, so results are unaffected. The Phase 4
  `porosity_bounded` case does change run behavior (it sets `0.065`); its
  paths/params were statically verified but the MOOSE solve was not executed here.
- Full compile/test execution not run in this worktree (no build tree; known
  macOS link caveat). The C++ change is string-literal-only and cannot alter
  behavior. Recommend a routine build + `test_ad_damaged_porosity.i` run on the
  cluster as a final gate.
- Reference table rows tagged `UNVERIFIED` (Heider & Markert 2017 page-level
  detail; the ~10–20% extreme secondary-porosity figure) must not be promoted to
  asserted fact without a resolving primary source.
