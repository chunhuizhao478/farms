# Implementation Plan: Bulk Modulus from the SPECTRAL Elasticity Tensor

## Overview
`NDSmallDeformationIsotropicElasticity` builds a damage- and strain-dependent
elastic tangent `C(d, ε)` under the **SPECTRAL** decomposition (the
decomposition used by every `2d_hydromech` production input). This plan adds a
new scalar material property `bulk_modulus_degraded` equal to the bulk modulus
**extracted from that elasticity tensor** by the volumetric contraction
`K = (1/9) I:C:I`, and surfaces it to Exodus output in the `2d_hydromech`
inputs. Scope is deliberately limited to the SPECTRAL decomposition — the
`NONE` and `VOLDEV` branches are out of scope.

## Definition of the quantity (math)

For a fourth-order elastic tangent `C = ∂σ/∂ε`, the bulk modulus is its
volumetric projection:

```
K_eff = (1/9) * I : C : I = (1/9) * Σ_{i,k} C_iikk ,   I = δ_ij (2nd-order identity)
```

We compute `K_eff` by **constructing the degraded SPECTRAL elasticity tensor
`C` and contracting it** — this is the literal "bulk modulus calculated from
the elasticity tensor." The SPECTRAL tangent (consistent with
`computeStressSpectralDecomposition`) is

```
C        = C_intact + (g − 1) · C_pos
C_intact = K·(I⊗I) + 2G·(I4_sym − (1/3) I⊗I)
C_pos    = λ·H(tr ε)·(I⊗I) + 2G·P⁺ ,    λ = K − 2G/LIBMESH_DIM
```

where `g = g(d)` is the degradation function, `P⁺ = ∂⟨ε⟩⁺/∂ε` is the
positive-projection tensor, `I⊗I = δ_ij δ_kl`, and `I4_sym` is the symmetric
fourth-order identity.

### Critical implementation constraint — do NOT reuse the existing Jacobian
The existing `computeJacobianSpectralDecomposition` builds `C` using
`RankFourTensor(initIdentity)` for the `I⊗I` (volumetric) terms. In this MOOSE
build `initIdentity` is the **diagonal-only** tensor `(*this)(i,i,i,i)=1`
(verified at `raccoon/moose/framework/include/utils/RankFourTensorImplementation.h:69-73`),
**not** `δ_ij δ_kl`. Contracting that tensor gives `I:(initIdentity):I = 3`
instead of `I:(I⊗I):I = 9`, so the extracted bulk modulus would be wrong
(`C_intact` would yield `(3K+4G)/9` instead of `K`). Therefore:

- Build the true `I⊗I` with **`RankTwoTensor::Identity().outerProduct(RankTwoTensor::Identity())`**
  (verified: `RankTwoTensor.h:1153`, `A.outerProduct(B)_ijkl = A_ij B_kl`).
- Do **not** call `computeJacobianSpectralDecomposition` for this purpose.
- Fixing the pre-existing `initIdentity` issue in the Jacobian itself is **out
  of scope** (it could change Newton convergence of production runs); it is
  recorded under Risks as a follow-up.

### Closed form (analytic cross-check, for test golds only)
Contracting the correctly-built `C` gives, using `I:(I⊗I):I = 9`,
`I:I4_sym:I = 3`, and `I:P⁺:I = Σ_a H(ε_a) = N⁺`:

```
K_eff = K + (g − 1) · ( λ·H_tr + (2G/9)·N⁺ )
```

- `H_tr = 1 if tr ε > 0 else 0`,
- `N⁺ = number of strictly positive principal strains (eig > 0, range 0..3)`.

Sanity checks:
- `d = 0` ⇒ `g = 1` ⇒ `K_eff = K` (intact), any strain.
- Full tension (`H_tr=1, N⁺=3`): `λ + 2G/3 = K` ⇒ `K_eff = g·K`.
- Full compression (`H_tr=0, N⁺=0`): `K_eff = K` (crack closed — volumetric
  stiffness preserved under compression).
- `K_eff > 0` strictly for `g ∈ [η, 1]` (bracket ≤ `K` ⇒ `K_eff ≥ g·K ≥ η·K`).

This formula is used only to derive expected test values; the implementation
computes `K_eff` from the assembled `C` tensor (so the contraction itself is
under test).

## Constraints
- **SPECTRAL only.** Do not implement `NONE`/`VOLDEV` variants of the property.
- **Build the tensor; do not reuse `computeJacobian*`.** Use the true `I⊗I`
  via `outerProduct` (see critical constraint above).
- Do not change the signatures of `computeStress` / `computeJacobian`, nor any
  existing material-property name or input-file parameter.
- The property must always be assigned (no declared-but-unassigned property);
  emit a one-time warning if `decomposition != SPECTRAL` since the value is
  defined for the spectral tangent.
- Follow existing style: hard-coded property name (cf. `effective_perm`,
  `crack_rotation`), `RankTwoTensor`/`RankFourTensor` MOOSE APIs, `LIBMESH_DIM`
  for `λ` (matches existing lines 427/519).
- Header member declaration order must match constructor initialization order
  (`-Wreorder`): append the new member **after** `_solid_bulk_compliance_damaged`.
- macOS link caveat (memory `toolchain-macos-tahoe-link-failure`): the library
  compiles even if the local `farms-opt` link fails; Phase-1 acceptance checks
  compilation, Phase-2 runs where the executable links.

---

## Phase 1: Compute the spectral bulk modulus in the material

### Goal
The material declares and populates `bulk_modulus_degraded = (1/9) I:C:I`,
computed by assembling the degraded SPECTRAL elasticity tensor `C` and
contracting it, at every quadrature point.

### Files to Create
- (none)

### Files to Modify
- `include/materials/phasefield_smalldeform_nonAD/NDSmallDeformationIsotropicElasticity.h`
  — declare the new property member and the helper method.
- `src/materials/phasefield_smalldeform_nonAD/NDSmallDeformationIsotropicElasticity.C`
  — declare the property in the initializer list, add the helper, assign the
  property in `computeStress`.

### Detailed Requirements

1. **Header — declare the helper** (private, near the other helpers, ~line 47):
   ```cpp
   // Bulk modulus extracted from the degraded SPECTRAL elastic tangent C(d, eps)
   // via the volumetric contraction K = (1/9) I:C:I. Builds C with the TRUE
   // I⊗I (outerProduct), NOT the diagonal-only RankFourTensor(initIdentity).
   Real computeSpectralBulkModulus(const RankTwoTensor & strain);
   ```

2. **Header — declare the new property member** as the **last** data member
   (immediately after `_solid_bulk_compliance_damaged`, line 178):
   ```cpp
   /// Bulk modulus K = (1/9) I:C:I of the degraded SPECTRAL elastic tangent.
   MaterialProperty<Real> & _bulk_modulus_degraded;
   ```

3. **Source — initializer list**: add as the **last** entry (after
   `_solid_bulk_compliance_damaged(...)`, line 254):
   ```cpp
   ,_bulk_modulus_degraded(declareProperty<Real>("bulk_modulus_degraded"))
   ```

4. **Source — helper implementation** (append after `getCharacteristicLength()`):
   ```cpp
   Real
   NDSmallDeformationIsotropicElasticity::computeSpectralBulkModulus(
       const RankTwoTensor & strain)
   {
     const Real g = _g[_qp];   // set by computeGDerivatives() before this call
     const Real K = _K[_qp];
     const Real G = _G[_qp];
     const Real lambda = K - 2.0 * G / LIBMESH_DIM;   // matches lines 427/519

     // TRUE second-order identity dyad I⊗I = δ_ij δ_kl. NOTE: do NOT use
     // RankFourTensor(initIdentity) — it is diagonal-only here, not I⊗I.
     const RankTwoTensor I2 = RankTwoTensor::Identity();
     const RankFourTensor IxI = I2.outerProduct(I2);
     const RankFourTensor I4_sym(RankFourTensor::initIdentitySymmetricFour);

     // Intact isotropic tangent: C_intact = K (I⊗I) + 2G (I4_sym − (1/3) I⊗I)
     const RankFourTensor C_intact = K * IxI + 2.0 * G * (I4_sym - IxI / 3.0);

     // Positive-projection part (spectral split), consistent with
     // computeStressSpectralDecomposition / computeJacobianSpectralDecomposition.
     RankTwoTensor eigvecs;
     std::vector<Real> eigvals(LIBMESH_DIM);
     const RankFourTensor P_pos =
         strain.positiveProjectionEigenDecomposition(eigvals, eigvecs);
     const Real H_tr = (strain.trace() > 0.0) ? 1.0 : 0.0;
     const RankFourTensor C_pos = lambda * H_tr * IxI + 2.0 * G * P_pos;

     // Degraded spectral elastic tangent.
     const RankFourTensor C = C_intact + (g - 1.0) * C_pos;

     // Bulk modulus = (1/9) I:C:I = (1/9) Σ_{i,k} C_iikk.
     Real ICI = 0.0;
     for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
       for (unsigned int k = 0; k < LIBMESH_DIM; ++k)
         ICI += C(i, i, k, k);
     return ICI / 9.0;
   }
   ```

5. **Source — assign the property in `computeStress`**: in `computeStress`
   (lines 369-387), immediately before `return stress;`, insert:
   ```cpp
   // Bulk modulus from the degraded SPECTRAL elastic tangent. g(d) is already
   // set by computeGDerivatives() above. Scope is the spectral model; warn once
   // if a different decomposition is active (the value assumes the spectral split).
   if (_decomposition != Decomposition::spectral)
     mooseDoOnce(mooseWarning(
         "bulk_modulus_degraded is defined for decomposition = SPECTRAL; the "
         "reported value assumes the spectral tangent."));
   _bulk_modulus_degraded[_qp] = computeSpectralBulkModulus(strain);
   ```

6. **(Optional, latent-bug cleanup — keep tightly scoped)**: the existing
   `solid_bulk_compliance_damaged` property is declared (line 254) but never
   assigned (its code is commented out at lines 655-661; no other file consumes
   it — verified by grep). If desired, assign it alongside the bulk modulus:
   ```cpp
   _solid_bulk_compliance_damaged[_qp] =
       1.0 / std::max(_bulk_modulus_degraded[_qp], 1e-30);
   ```
   and delete the dead comment block at lines 655-661. This is **optional** and
   may be deferred to keep this change focused on the bulk-modulus output.

### Interfaces
- New material property `bulk_modulus_degraded` (`MaterialProperty<Real>`),
  readable by `MaterialRealAux`, `ElementAverageValue` (via aux),
  `ElementIntegralMaterialProperty`, or the material `output_properties`.
- New private method `Real computeSpectralBulkModulus(const RankTwoTensor &)`.

### Edge Cases to Handle
- **Undamaged (`d=0`, `g=1`)** → `C = C_intact` → `K_eff = K`, any strain.
- **Full compression** (all eigenvalues ≤ 0) → `H_tr=0`, `P⁺=0` → `K_eff = K`
  (crack closure preserves volumetric stiffness).
- **2D plane strain out-of-plane eigenvalue**: `ComputeSmallStrain` yields
  `ε_zz = 0`; `0` is **not** counted in `N⁺` (strict `>0`), so equibiaxial
  tension gives `N⁺=2` (intermediate), not `3`. Documented in test golds.
- **Non-spectral decomposition**: property still assigned (spectral
  interpretation) + one-time warning; no unassigned-property bug.

### Acceptance Criteria
- [ ] Header + source compile cleanly (`make -j`; object builds even if the
      local `farms-opt` link fails per the macOS caveat).
- [ ] No `-Wreorder` warning (member declared/initialized last).
- [ ] Implementation builds `IxI` via `outerProduct` and contracts `C(i,i,k,k)`;
      it does **not** call `computeJacobianSpectralDecomposition`.
- [ ] `bulk_modulus_degraded == K` for `d=0` and for SPECTRAL compression
      (checked in Phase 2).

### Dependencies
- Depends on: nothing.
- Required by: Phase 2 (test), Phase 3 (output wiring).

---

## Phase 2: Single-element regression test (SPECTRAL)

### Goal
A deterministic, analytically-checked test proves the value extracted from the
spectral tensor matches the closed form across the three regimes
(undamaged → `K`, damaged compression → `K`, damaged tension → degraded).

### Files to Create
- `test/tests/materials/permeability_normal_strain/bulk_modulus_spectral.i`
  — single QUAD4, function-driven `d`/`disp`, reads `bulk_modulus_degraded`
  into an aux var + `ElementAverageValue` postprocessor; CSV output. Mirror the
  structure of `isotropic_axis.i` (same `[Mesh]`, trivial `[u]` diffusion solve,
  `FunctionAux` for `d`/`disp_x`/`disp_y`).

### Files to Modify
- `test/tests/materials/permeability_normal_strain/tests` — add three CSVDiff
  blocks.
- `test/tests/materials/permeability_normal_strain/EXPECTED_VALUES.md` — add a
  "Degraded bulk modulus (SPECTRAL)" section with the derivations.

### Detailed Requirements
1. **Input** `bulk_modulus_spectral.i`:
   - Reuse `[Mesh]`/`[GlobalParams]`/`[Variables u]`/`[Kernels diff]`/`[BCs]`
     from `isotropic_axis.i`.
   - `[Materials]`: `ComputeIsotropicElasticityTensor`
     (`youngs_modulus=50e9 poissons_ratio=0.3`), `ComputeSmallStrain`,
     `GenericConstantMaterial` `K G = '4.1666666667e10 1.9230769231e10'`,
     `NDSmallDeformationIsotropicElasticity` with `decomposition=SPECTRAL
     model_type=AT1 eta=1e-6 porous_flow_coupling=false`, and
     `NDComputeSmallDeformationStress`. Name the `[elasticity]` block so it can
     be referenced from `cli_args`.
   - One input driven by `cli_args` (no duplication); tests set `d`/`disp`
     functions and `file_base`.
   - `[AuxVariables] Kd` (`MONOMIAL CONSTANT`); `[AuxKernels]`
     `MaterialRealAux property=bulk_modulus_degraded variable=Kd
     execute_on='TIMESTEP_END'`, plus `FunctionAux` `d_aux`, `disp_x_aux`,
     `disp_y_aux`.
   - `[Postprocessors] Kd_avg type=ElementAverageValue variable=Kd
     execute_on='INITIAL TIMESTEP_END'` (unit element ⇒ average == value).
   - `[Outputs] csv = true`.

2. **`tests` entries** (`K=4.1666666667e10`, `G=1.9230769231e10`, `η=1e-6`,
   AT1 `g=(1−d)²(1−η)+η`, `λ=K−2G/3=2.8846153846e10`):
   ```
   [bulk_modulus_spectral_undamaged]
     # d=0 -> g=1 -> K_eff = K regardless of strain (tension imposed).
     type = 'CSVDiff'
     input = 'bulk_modulus_spectral.i'
     cli_args = 'AuxKernels/d_aux/function=0
                 AuxKernels/disp_x_aux/function=1e-3*x
                 AuxKernels/disp_y_aux/function=1e-3*y
                 Outputs/file_base=bulk_modulus_spectral_undamaged'
     csvdiff = 'bulk_modulus_spectral_undamaged.csv'
     rel_err = 1e-8
     requirement = 'The system shall extract K_eff = K from the SPECTRAL elastic '
                   'tensor in the undamaged limit (g=1).'
   []

   [bulk_modulus_spectral_compression]
     # d=0.7, equibiaxial compression -> H_tr=0, N+=0 -> K_eff = K (closure).
     type = 'CSVDiff'
     input = 'bulk_modulus_spectral.i'
     cli_args = 'AuxKernels/d_aux/function=0.7
                 AuxKernels/disp_x_aux/function=-1e-3*x
                 AuxKernels/disp_y_aux/function=-1e-3*y
                 Outputs/file_base=bulk_modulus_spectral_compression'
     csvdiff = 'bulk_modulus_spectral_compression.csv'
     rel_err = 1e-8
     requirement = 'The system shall extract the full bulk modulus K_eff = K from '
                   'the SPECTRAL tensor under pure compression at d=0.7 (crack closure).'
   []

   [bulk_modulus_spectral_tension]
     # d=0.7, equibiaxial tension eps=diag(1e-3,1e-3,0) -> H_tr=1, N+=2.
     # K_eff = K + (g-1)*(lambda + (2G/9)*2) ~= 7.638923e9.
     type = 'CSVDiff'
     input = 'bulk_modulus_spectral.i'
     cli_args = 'AuxKernels/d_aux/function=0.7
                 AuxKernels/disp_x_aux/function=1e-3*x
                 AuxKernels/disp_y_aux/function=1e-3*y
                 Outputs/file_base=bulk_modulus_spectral_tension'
     csvdiff = 'bulk_modulus_spectral_tension.csv'
     rel_err = 1e-8
     requirement = 'The system shall extract the degraded bulk modulus from the '
                   'SPECTRAL tensor under tension at d=0.7, matching the analytic '
                   'K + (g-1)(lambda + (2G/9)N+) with N+=2 (out-of-plane eig=0).'
   []
   ```
   (Confirm the exact `cli_args` paths against the block/object names used.)

3. **Expected golds** (hand-computed; place gold CSVs under `gold/`):
   - undamaged: `K_eff = K = 4.1666666667e10`.
   - compression: `K_eff = K = 4.1666666667e10`.
   - tension: `g = 0.09000091`, `bracket = λ + (2G/9)·2 = 3.7393162393e10`,
     `K_eff = K + (g−1)·bracket = 7.638923e9`.

4. **EXPECTED_VALUES.md**: append the SPECTRAL bulk-modulus derivation
   (`K_eff = (1/9) I:C:I`, the closed form, and the three regime values),
   matching the existing documentation style.

### Edge Cases to Handle
- `porous_flow_coupling=false` must still produce `bulk_modulus_degraded`
  (decoupled from permeability).
- The tension case is the discriminating test: it fails if `IxI` is built from
  the diagonal-only `initIdentity` (would shift the volumetric contraction) or
  if `P⁺` is mis-contracted.

### Acceptance Criteria
- [ ] All three CSVDiff cases pass within `rel_err=1e-8`.
- [ ] Compression and undamaged both give `K` (asymmetry vs. the tension case
      `7.638923e9` at the same `d=0.7` is the key physics check).
- [ ] Existing `permeability_normal_strain` tests still pass unchanged.

### Dependencies
- Depends on: Phase 1.
- Required by: nothing.

---

## Phase 3: Surface K_eff in the 2d_hydromech production inputs

### Goal
`bulk_modulus_degraded` is written to Exodus in the `2d_hydromech` inputs
(all of which use `decomposition = SPECTRAL`).

### Files to Modify (canonical, do first)
- `pulsepower/cmame_revision/2d_hydromech/permeability_formula/elasticity_E1d25.i`

Minimal `output_properties` wiring (the `[elasticity]` block already uses it,
line 670):
1. Extend `output_properties = 'elastic_strain psie_active'` to
   `output_properties = 'elastic_strain psie_active bulk_modulus_degraded'`
   (`outputs = exodus` already present; the field is auto-created as
   `CONSTANT MONOMIAL` named `bulk_modulus_degraded`).
2. Append `bulk_modulus_degraded` to the `[Outputs]/[exodus]` `show = '...'`
   list (line 909 — explicit list, so the field is hidden unless named).

> Alternative (explicit aux field, mirroring `biot_modulus_aux`): add an
> `[AuxVariables] bulk_modulus_degraded_aux` (`MONOMIAL CONSTANT`), a
> `MaterialRealAux` (`property=bulk_modulus_degraded execute_on='TIMESTEP_END'`),
> and add it to `show`. Pick one; default is the `output_properties` approach.

### Files to Modify (mechanical replication — same two-line pattern)
Apply the identical edit to the other `NDSmallDeformationIsotropicElasticity`
inputs under `2d_hydromech`:
- `permeability_formula/pulse_duration_magnitude/elasticity_E1d25_pulse4em5.i`
- `permeability_formula/pulse_duration_magnitude/elasticity_E1d25_pulse2em5.i`
- `permeability_formula/pulse_duration_magnitude/elasticity_E1d25.i`
- `permeability_formula/pulse_duration_magnitude/elasticity_E2d50.i`
- `permeability_formula/pulse_duration_magnitude/elasticity_E5d00.i`
- `permeability_formula/undrained/em_0p005/elasticity_E1d25.i`
- `permeability_formula/undrained/em_0p010/elasticity_E1d25.i`
- `permeability_formula/undrained/em_0p020/elasticity_E1d25.i`
- `permeability_formula/undrained/em_0p040/elasticity_E1d25.i`
- `permeability_formula/undrained/em_0p080/elasticity_E1d25.i`
- `permeability_formula/drained/em_0p005/elasticity_E1d25.i`
- `permeability_formula/drained/em_0p010/elasticity_E1d25.i`
- `permeability_formula/drained/em_0p020/elasticity_E1d25.i`
- `permeability_formula/drained/em_0p080/elasticity_E1d25.i`
- `permeability_formula/order_test/em_0p005/elasticity_E1d25.i`
- `permeability_formula/viscosity/mu_1em3/elasticity_E1d25.i`
- `permeability_formula/viscosity/mu_1em4/elasticity_E1d25.i`
- `permeability_formula/oil_water_mixture/elasticity_E1d25.i`
- `benchmark/elasticity_E1d25.i`

> Scope note: mechanical and low-risk but ~20 files. Do the canonical file +
> one run first, confirm the field appears, then replicate. Do not add the
> property to inputs that do not instantiate the material (e.g. `fracture_*.i`,
> `static_solve.i`).

### Edge Cases to Handle
- Each input must still parse (`--check-input`). The new field is diagnostic
  only and must not alter physics, residuals, or energy postprocessors.

### Acceptance Criteria
- [ ] `permeability_formula/elasticity_E1d25.i` parses (`--check-input` clean).
- [ ] A short run writes a non-trivial `bulk_modulus_degraded` field
      (≈ `K` in undamaged/compressed regions, decreasing in damaged tensile zones).
- [ ] No change to energy-balance postprocessor output vs. a pre-change run.

### Dependencies
- Depends on: Phase 1. Phase 2 recommended first.
- Required by: nothing.

---

## Testing Strategy
- **Phase 1**: compile check; reasoning against the closed form.
- **Phase 2**: three CSVDiff regressions with analytic golds (undamaged → `K`,
  compression → `K`, tension → `7.638923e9`), all exercising the actual tensor
  contraction.
- **Regression**: rerun the full `permeability_normal_strain` suite (no
  `effective_perm` behavior change).
- **Validation**: in a production run, confirm `bulk_modulus_degraded → K`
  where `d→0`/compressed and decreases in damaged tensile zones.

## Risk Assessment
- **`initIdentity` is diagonal-only** (`RankFourTensorImplementation.h:69`).
  *Mitigation*: build `I⊗I` with `outerProduct`; never reuse the existing
  Jacobian. This is the single most important constraint — the Phase-2 tension
  test is designed to catch a regression to the diagonal tensor.
- **Pre-existing Jacobian inconsistency**: `computeJacobianSpectralDecomposition`
  uses the diagonal `initIdentity` for its `I⊗I` terms, so the analytic Jacobian
  likely does not equal the true tangent of the spectral stress. **Out of scope**
  here (could affect Newton convergence of production runs); recommend a
  follow-up ticket. Do not let `K_eff` depend on it.
- **`P⁺` contraction**: the implementation relies on `I:C:I` via the assembled
  tensor, so correctness does not depend on the analytic `I:P⁺:I = N⁺` identity
  — that identity is only used to hand-compute the gold. If the tension gold and
  the run disagree, suspect the strain the material sees (plane-strain `ε_zz`)
  rather than the contraction.
- **Member ordering**: append the new member/initializer last (`-Wreorder`).
- **macOS local link failure** (memory `toolchain-macos-tahoe-link-failure`):
  Phase-1 acceptance is compilation; run Phase-2 where `farms-opt` links.
- **Phase-3 breadth**: ~20 near-identical inputs; a `show`-list typo silently
  hides the field — verify on the canonical file first.
