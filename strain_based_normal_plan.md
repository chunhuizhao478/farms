# Implementation Plan: Strain-Based Crack Normal (Liu et al. 2024, CMAME Eqs. 29–30)

## Overview
Rewrite **only the crack-normal calculation** of the existing `normal_strain`
permeability model so the unit crack normal is the eigenvector of the **maximum
principal strain**, `n_F = e₁` (Liu et al. 2024, *Comput. Methods Appl. Mech.
Engrg.* **429**:117165, Eqs. 29–30). We do **not** touch the permeability assembly
(Eq. 28 is explicitly out of scope): the Heider-2021 aperture `w_c = h_c·|1+ε_nn|`,
the `d^b` weighting, the tangential projector `(I − n⊗n)`, the Heaviside gate, the
residual aperture, and the whole `regularize_crack_normal` damage-gradient path all
stay byte-for-byte unchanged. The strain-based normal is then A/B compared against
the existing `regularize_crack_normal` (damage-gradient) normal under
`pulsepower/cmame_revision/2d_hydromech`, isolating the effect of the **normal
definition** alone.

## Scope correction vs. the previous draft
- **Dropped:** new `permeability_model = strain_based`, Eq. 28 reassembly, Eq. 31/32
  aperture replacement, and the `aperture_strain_source` / `clamp_aperture_strain`
  parameters. None of that is needed.
- **Kept / refined:** the strain-based normal `n_F = e₁` (Eqs. 29–30). This is the
  existing `crack_normal_source = principal_strain` option, made faithful to Eq. 29
  by deriving `e₁` from the model's **mechanical/total strain** (`_total_strain`)
  rather than from `_elastic_strain` via `_crack_rotation`.

## The Math (what changes)
- **Eq. 29 (principal decomposition):** `ε = Σ_i ε_i e_i`, `ε₁ ≥ ε₂ ≥ ε₃`, where `ε`
  is the (total/mechanical) linearized strain.
- **Eq. 30 (crack normal):** `n_F = e₁` — eigenvector of the **largest** principal
  strain. `RankTwoTensor::symmetricEigenvaluesEigenvectors(eigval, eigvec)` returns
  eigenvalues in **ascending** order, so `ε₁ = eigval[2]` and `e₁ = eigvec.column(2)`;
  eigenvectors are unit-norm by construction.

Everything downstream is unchanged and reproduced here only for context:
`ε_nn = n_F·ε·n_F`, `w_c = h_c·|1+ε_nn|`, `χ_d = H(d − d_thr)`,
`w_h = max(f_c·w_c·χ_d, f_c·w_r·χ_d)`, `k_w = w_h²/12`,
`K = k₀·I + d^b·k_w·(I − n_F⊗n_F)`.

## Why this is a real (small) code change, not a no-op
`crack_normal_source = principal_strain` already returns `e₁`, **but** from
`_crack_rotation`, which `computeCrackStrainAndOrientation` (`.C:664–691`) builds by
eigendecomposing `_elastic_strain`. The aperture term `ε_nn` instead reads
`_total_strain` (= `mechanical_strain`, `.C:810–813`). In pure elasticity these
strains are identical, but they diverge if a plastic strain is ever introduced, and
Eq. 29 is defined on the total strain. The rewrite computes `n_F` from the **same**
`_total_strain` used for `ε_nn`, making the normal calculation self-consistent and
literally Eq. 29–30. Because the unit tests are pure-elastic, **existing gold files
do not change** — this is verified, not assumed (see Phase 2 regression criterion).

## Constraints
- **Do not change the permeability assembly or any other branch.** The
  `none / exponential / darcy_poiseuille` branches, the `damage_gradient` normal,
  `regularize_crack_normal`, the aperture/threshold/residual logic, and `d^b`
  weighting must be untouched so the A/B baseline is preserved.
- **Reuse existing members.** `_total_strain` (already bound under `normal_strain`,
  `.C:314–315`) and the existing eigen-API. No new input parameters unless noted.
- **Keep `_crack_rotation` intact.** It is still produced by
  `computeCrackStrainAndOrientation` and still consumed by the legacy
  exponential/Darcy `rotate(R)` calls; only the `principal_strain` *normal* stops
  reading it.
- **Execution path.** `updatePermeabilityForCracking()` runs only under
  `decomposition = SPECTRAL` (unchanged, pre-existing limitation shared by the whole
  `normal_strain` model). The production input already uses `SPECTRAL`.
- **Conventions.** Non-AD `RankTwoTensor`/`RealVectorValue`; `paramError` for
  validation; comments cite "Liu et al. 2024 CMAME Eq. (29)/(30)".

---

## Phase 1: Rewrite the `principal_strain` crack normal to be strain-based (Eqs. 29–30)

### Goal
After this phase, `crack_normal_source = principal_strain` computes
`n_F = e₁` from the model's mechanical/total strain `_total_strain` (Eqs. 29–30),
self-consistent with `ε_nn`, while every other code path and all 8 existing gold
files are unchanged.

### Files to Modify
- `include/materials/phasefield_smalldeform_nonAD/NDSmallDeformationIsotropicElasticity.h`
  — declare a small private helper.
- `src/materials/phasefield_smalldeform_nonAD/NDSmallDeformationIsotropicElasticity.C`
  — add the helper body; replace the 3-line `principal_strain` normal read; refresh
  the `crack_normal_source` parameter doc-string to cite Eqs. 29–30.

### Detailed Requirements

**1.1 — Add a private helper (header, near `getCharacteristicLength` decl, ~line 47).**
```cpp
// Liu et al. 2024 CMAME eqs. (29)-(30): return the unit eigenvector e_1 of the
// MAXIMUM principal strain of `strain` (the strain-based crack normal n_F).
RealVectorValue maxPrincipalStrainDirection(const RankTwoTensor & strain) const;
```

**1.2 — Implement the helper (`.C`, place next to `spectralDecomposition`,
after ~line 589).**
```cpp
RealVectorValue
NDSmallDeformationIsotropicElasticity::maxPrincipalStrainDirection(
    const RankTwoTensor & strain) const
{
  // Eq. (29): eps = sum_i eps_i e_i, eigenvalues ascending.
  std::vector<Real> eigval(3, 0.0);
  RankTwoTensor eigvec;
  strain.symmetricEigenvaluesEigenvectors(eigval, eigvec);
  // Eq. (30): n_F = e_1 = eigenvector of the largest principal strain.
  // Ascending order => column(2) is the most-tensile direction; unit-norm.
  return eigvec.column(2);
}
```
(`symmetricEigenvaluesEigenvectors` is const-callable — it is already invoked on the
const `_elastic_strain[_qp]` at `.C:679` — so this method is `const`.)

**1.3 — Replace the `principal_strain` normal read (`.C:769–780`).**

Replace the existing block:
```cpp
    else // principal_strain fallback
    {
      // _crack_rotation[_qp] column 0 is the most-tensile eigenvector, ...
      n_d(0) = R(0, 0);
      n_d(1) = R(1, 0);
      n_d(2) = R(2, 0);
      have_normal = true;
    }
```
with:
```cpp
    else // principal_strain: strain-based crack normal, Liu 2024 eqs. (29)-(30)
    {
      // Eq. (29)-(30): n_F = e_1 = eigenvector of the maximum principal strain.
      // Derive it from the model's mechanical strain `_total_strain` (the same
      // strain used for eps_nn below), NOT from _crack_rotation/_elastic_strain,
      // so the normal and the normal-strain aperture are computed from one
      // self-consistent strain tensor (identical in pure elasticity; the correct
      // total strain if plasticity is later attached). Eigenvectors are unit-norm,
      // so n_d is a valid unit normal everywhere (no |grad d| division, no
      // regularization needed -- contrast crack_normal_source = damage_gradient).
      n_d = maxPrincipalStrainDirection((*_total_strain)[_qp]);
      have_normal = true;
    }
```
No other lines in `updatePermeabilityForCracking()` change. `R = _crack_rotation[_qp]`
(`.C:702`) remains for the exponential/Darcy `rotate(R)` branches.

**1.4 — Refresh the `crack_normal_source` parameter doc-string (`.C:79–84`).**

Update the description so `principal_strain` cites the equations and the strain
source (no enum values added/removed):
```cpp
params.addParam<MooseEnum>(
    "crack_normal_source",
    MooseEnum("damage_gradient principal_strain", "damage_gradient"),
    "Source of the unit crack normal n_F. 'damage_gradient' uses "
    "grad(d)/|grad(d)| (Heider 2021 eq. 46), optionally regularized via "
    "regularize_crack_normal. 'principal_strain' uses the strain-based normal "
    "n_F = e_1, the eigenvector of the maximum principal strain of the model's "
    "mechanical strain (Liu et al. 2024 CMAME eqs. 29-30); it is unit-norm "
    "everywhere and needs no regularization.");
```

### Interfaces
- New private `RealVectorValue maxPrincipalStrainDirection(const RankTwoTensor &) const`.
- No new parameters, no enum changes, no member additions beyond the helper.
- Consumes existing `_total_strain` and the existing eigen-API.

### Edge Cases to Handle
- **Zero strain** (`ε = 0`): eigenvectors are arbitrary but `ε_nn = 0` ⇒
  `w_c = h_c` and the projector still well-defined; result is finite. This matches
  prior `principal_strain` behavior (it also returned an arbitrary unit vector for a
  degenerate tensor). No regression.
- **Degenerate `ε₁ = ε₂`**: `e₁` is non-unique but still unit; pick non-degenerate
  strains in the discriminating test (Phase 2).
- **`_total_strain == nullptr`**: cannot occur — it is bound whenever
  `permeability_model = normal_strain` (`.C:314–315`), the only model that reaches
  this branch. The constructor already `paramError`s if `mechanical_strain` is
  absent (`.C:308–313`).
- **Fully-damaged core** (`d → 1`, `grad d → 0`): unlike the damage-gradient normal,
  the strain-based normal stays well-defined (`e₁` of the strain), so the anisotropic
  projector is retained at the core — this is precisely the effect the A/B study
  isolates against `regularize_crack_normal` (which drives `n_d → 0`, isotropic core).

### Acceptance Criteria
- [ ] Compiles (`make -j`) with no new warnings.
- [ ] All 8 existing tests in `test/tests/materials/permeability_normal_strain/`
      pass with **unchanged** gold files (pure-elastic ⇒ total = elastic strain).
- [ ] In a case where `grad(d)` and `e₁` point in **different** directions,
      `crack_normal_source = principal_strain` yields a normal along `e₁` (Phase 2
      test `strain_based_normal.i`), not along `grad(d)`.

### Dependencies
- Depends on: nothing.
- Required by: Phase 2, Phase 3.

---

## Phase 2: Discriminating unit test + regression confirmation

### Goal
Pin that the rewritten `principal_strain` normal follows the **strain** (`e₁`), not
the damage gradient, with a single-QUAD4 test where the two normals differ; confirm
the existing gold files are untouched.

### Files to Create
- `test/tests/materials/permeability_normal_strain/strain_based_normal.i`

### Files to Modify
- `test/tests/materials/permeability_normal_strain/tests` — add one `Exodiff` block.
- `test/tests/materials/permeability_normal_strain/EXPECTED_VALUES.md` — append the
  hand calc below.

### Test design (`strain_based_normal.i` — copy `principal_strain_fallback.i`, edit aux + material)
Construct a case where the damage gradient and the max-principal-strain direction are
**orthogonal**, so a damage-based normal and a strain-based normal give different
effective-permeability tensors:
- Damage: `d(x,y) = 0.5 + 0.4*x` ⇒ `grad(d) = (0.4, 0, 0)` ⇒ damage-gradient normal
  would be `e_x`. (`d_qp = 0.5844, 0.8156` per the 2×2 Gauss points.)
- Strain: `disp_x = 0`, `disp_y = 1e-3*y` ⇒ `ε = diag(0, 1e-3, 0)` ⇒ `ε₁ = 1e-3`,
  `e₁ = e_y`.
- `[elasticity]` keys: `permeability_model = normal_strain`,
  `crack_normal_source = principal_strain`, `intrinsic_permeability = 5e-19`,
  `perm_exponent = 2`, `characteristic_length_type = constant`,
  `characteristic_length_value = 1e-3`, `permeability_anisotropic = true`,
  `damage_threshold_for_permeability = 0.5`, `correction_factor_fc = 1.0`.

Hand calc (strain-based normal `n_F = e_y`):
- `ε_nn = n_F·ε·n_F = ε_yy = 1e-3`; `w_c = 1e-3·1.001 = 1.001e-3`;
  `k_w = (1.001e-3)²/12 = 8.350e-8`.
- Projector `I − e_y⊗e_y = diag(1, 0, 1)`.
- `⟨d²⟩ = ½(0.5844² + 0.8156²) = 0.5033`; `α = ⟨d²⟩·k_w = 0.5033·8.350e-8 ≈ 4.202e-8`.
- **Expected:** `K_yy = k₀ = 5e-19` (normal direction blocked), `K_xx = K_zz = α ≈
  4.202838e-8`, `K_xy = 0`.
- **Discriminator:** a damage-gradient normal (`e_x`) would instead give
  `K_xx = k₀`, `K_yy = K_zz = α`. The swap of which diagonal entry collapses to `k₀`
  is what proves the normal is strain-based, not damage-based.

`tests` block (follow existing style; `design`/`issues` inherited):
```
[strain_based_normal]
  type = 'Exodiff'
  input = 'strain_based_normal.i'
  exodiff = 'strain_based_normal_out.e'
  requirement = 'With crack_normal_source = principal_strain, the system shall '
                'compute the crack normal as n_F = e_1, the eigenvector of the '
                'maximum principal strain (Liu 2024 CMAME eqs. 29-30), taken from '
                'the mechanical strain, so that when grad(d) || e_x but e_1 = e_y '
                'the tangential projector blocks the y-direction (K_yy = k0) and '
                'enhances x,z -- distinct from the damage-gradient normal.'
[]
```

### Gold regeneration
```bash
cd test/tests/materials/permeability_normal_strain
../../../../farms-opt -i strain_based_normal.i     # or the project run wrapper
cp strain_based_normal_out.e gold/
```
If the executable cannot be linked in the implementing environment, commit the input
+ `EXPECTED_VALUES.md` entry and flag the gold `.e` as **pending regeneration** (same
precedent as `regularized_normal_core`, EXPECTED_VALUES.md Test 8).

### Acceptance Criteria
- [ ] `strain_based_normal` passes Exodiff against gold, matching the hand calc
      (`K_yy = k₀`, `K_xx = K_zz ≈ 4.20e-8`).
- [ ] The 8 pre-existing tests still pass with unchanged gold (regression guard for
      "total = elastic strain in pure elasticity").
- [ ] `EXPECTED_VALUES.md` documents the discriminator vs. the damage-gradient normal.

### Dependencies
- Depends on: Phase 1.
- Required by: nothing.

---

## Phase 3: A/B experiment in `pulsepower/cmame_revision/2d_hydromech`

### Goal
Run the production case with the **only** difference being the crack-normal
definition — strain-based (`principal_strain`, Eqs. 29–30) vs. the existing
regularized damage-gradient normal — so the effect on fracture permeability and the
resulting hydro-mechanical response can be compared.

### Files to Create
- `permeability_formula/strain_based_normal/elasticity_E1d25.i` — copy of
  `permeability_formula/elasticity_E1d25.i` with the `[elasticity]` block's normal
  switched (diff below).
- `permeability_formula/strain_based_normal/fracture_E1d25.i`,
  `permeability_formula/strain_based_normal/static_solve.i` — verbatim copies of the
  siblings (the staggered sub-app and static initializer do not reference the
  permeability normal); keep relative MultiApp/Transfer paths consistent.
- `permeability_formula/strain_based_normal/submit_elasticity.sbatch`,
  `submit_static.sbatch`, `sync_results.sh` — copies adjusted for the new path.
- (Recommended baseline sibling) `permeability_formula/regularized_normal/elasticity_E1d25.i`
  — the existing input with `regularize_crack_normal = true` made explicit, so the
  comparison is "strain-based normal vs. regularized damage-gradient normal".

### `[elasticity]` block diff — strain-based variant
From the current production block (`permeability_formula/elasticity_E1d25.i:676–685`),
change only the normal source:
```
##----- normal_strain permeability, STRAIN-BASED normal (Liu 2024 eqs. 29-30) -----##
permeability_model = normal_strain
intrinsic_permeability = ${intrinsic_permeability}
perm_exponent = ${perm_exponent}
crack_normal_source = principal_strain    # n_F = e_1 (max principal strain), eqs. 29-30
characteristic_length_type = element_size # h_c = element size (paper default)
element_size_variable = mesh_size         # reuse existing ElementLengthAux aux var
permeability_anisotropic = true           # K_frac = (w_c^2/12)(I - n_F (x) n_F)
damage_threshold_for_permeability = 0.5
correction_factor_fc = 1.0
##---------------------------------------------------------------------------------##
```
Remove `regularize_crack_normal` / `crack_normal_regularization` / `damage_gradient`
keys from this variant. Everything else (mesh, `[strain] ComputeSmallStrain`,
`mesh_size` aux + `ElementLengthAux`, porous-flow kernels, history-energy material)
stays identical to the baseline.

### `[elasticity]` block — regularized-normal baseline (A side)
Existing production block plus the explicit regularization keys:
```
crack_normal_source = damage_gradient
regularize_crack_normal = true
crack_normal_regularization = 1e-8      # match the value used in the 3d cases
```

### Acceptance Criteria
- [ ] Both variants parse (`--check-input`) with no unused-parameter warnings; a
      `diff` of the two `[elasticity]` blocks shows only the normal-source keys differ.
- [ ] A short run (a few steps) of each writes `effective_perm` to Exodus, enabling
      side-by-side comparison of the fracture-permeability field along the crack.
- [ ] Qualitative check: at the fully-damaged crack core the strain-based variant
      retains an anisotropic projector (normal-direction conductivity suppressed),
      whereas the `regularize_crack_normal` baseline trends isotropic — the expected
      distinguishing effect.

### Dependencies
- Depends on: Phase 1.
- Required by: nothing.

---

## Testing Strategy
- **Unit (Phase 2):** one discriminating single-QUAD4 Exodiff where `grad(d) ⟂ e₁`,
  validated against a closed-form hand calc; the 8 existing gold files act as the
  regression guard that nothing else moved.
- **A/B (Phase 3):** field comparison of `effective_perm` between the strain-based
  normal and the regularized damage-gradient normal on the production mesh.
- **Optional sign check:** a rotated-principal-direction case (reuse the
  `anisotropic_rotated` CSV pattern) if a sign-sensitive projector test is wanted for
  the strain-based normal.

## Risk Assessment
- **Eigenvalue ordering.** Ascending ⇒ `e₁ = eigvec.column(2)`. Using `column(0)`
  would pick the most-compressive direction. *Detection:* `strain_based_normal.i`
  asserts `K_yy = k₀` (not `K_xx`).
- **Hidden gold drift.** If any reviewer suspects the rewrite changes elastic-case
  results, re-run the 8 existing tests; total = elastic strain in pure elasticity
  guarantees identical output. *Detection:* the regression criterion in Phase 1/2.
- **`_total_strain` availability.** Bound only under `normal_strain`; the
  `principal_strain` branch is only reachable under `normal_strain`, so it is always
  bound. *Detection:* immediate crash if mis-wired; the constructor's
  `mechanical_strain` `paramError` fires first when `ComputeSmallStrain` is absent.
- **Const-correctness of the helper.** `symmetricEigenvaluesEigenvectors` must be
  const-callable; it already is (invoked on `const _elastic_strain[_qp]`).
  *Detection:* compile error if not.
- **Decomposition requirement.** The normal is only updated under
  `decomposition = SPECTRAL` (pre-existing). The production input uses `SPECTRAL`;
  note it in the experiment README. *Detection:* `effective_perm` stays at its
  initial value if a non-SPECTRAL decomposition is used.
