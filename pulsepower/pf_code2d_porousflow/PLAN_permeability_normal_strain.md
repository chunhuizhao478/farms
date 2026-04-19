# Implementation Plan: Normal-Strain-Driven Permeability Enhancement

## Overview

Replace the simplified damage-based crack aperture `w = d * wc` used in the current
Darcy–Poiseuille permeability model with a kinematically consistent aperture derived
from the **normal strain to the fracture plane**, following the Heider (2021) review
of phase-field hydraulic fracturing (eqs. 46–48). The crack normal is computed from
the damage gradient `n_d = ∇d / |∇d|`; the aperture is `w_c = h_c · |1 + n_d · ε · n_d|`
(Heider eq. 47), where `h_c` is a characteristic length (element size `h`,
regularization length `l`, or a user constant). Fracture permeability is anisotropic,
enhancing flow in the crack-parallel plane via the tangential projector
`I − n_d ⊗ n_d`, and is blended with the matrix permeability through a damage-weighted
form `K = K_poro + (d^S)^b · K_frac` with a Heaviside gate `χ_d = H(d − d_threshold)`
(default `d_threshold = 0.5`, per the paper). The change is implemented as an
additional permeability model option in `NDSmallDeformationIsotropicElasticity`, with
a corresponding update of the two input files cited by the user.

## Constraints

### Interface constraints (cannot change)
- Public MOOSE material property **`effective_perm`** (`MaterialProperty<RealTensorValue>`)
  is the sole output consumed downstream by
  `ElkPorousFlowPermeabilityDamaged::computeQpProperties()`
  (`src/materials/porousflowmatprops/ElkPorousFlowPermeabilityDamaged.C:22`).
  Its type, name, and semantics (absolute permeability tensor in global frame,
  SI units m²) must not change.
- `ElkPorousFlowPermeabilityDamaged` and every kernel tagged with `use_damaged_biot = true`
  in `elasticity_E1d25.i` (poro_x, poro_y, mass0, flow_fluid_driving_energy) already
  read material properties (`biot_coefficient_damaged`, `PorousFlow_porosity_qp_damaged`,
  `PorousFlow_constant_biot_modulus_qp`) — none of their interfaces may change.
- `computeCrackStrainAndOrientation(RealVectorValue &)` and
  `updatePermeabilityForCracking()` are private virtual methods of
  `NDSmallDeformationIsotropicElasticity`
  (`include/materials/phasefield_smalldeform_nonAD/NDSmallDeformationIsotropicElasticity.h:37–38`).
  Their contract with `computeStressSpectralDecomposition`
  (`src/materials/phasefield_smalldeform_nonAD/NDSmallDeformationIsotropicElasticity.C:216–253`)
  is the only call site and may be modified internally, but they must remain
  side-effect compatible (write to `_crack_rotation` and `_effective_perm`).

### Dependency constraints
- MOOSE `RankTwoTensor::symmetricEigenvaluesEigenvectors(eigval, eigvec)` — already used.
- MOOSE `RankTwoTensor::rotate(R)` — already used; performs `R · A · Rᵀ` in place.
- `_elastic_strain[_qp]` (inherited from `NDSmallDeformationElasticityModel`) — already used.
- No new MOOSE core dependencies.

### Convention constraints
- Follow the existing pattern used for `exponential_permeability_model` and
  `darcy_poiseuille_permeability_model`: add a bool flag plus its own companion
  Real parameters, validate mutual exclusivity in the constructor, and branch
  inside `updatePermeabilityForCracking()`.
- Keep the MooseObject class registration macro and file organization unchanged
  (header in `include/materials/phasefield_smalldeform_nonAD/`, implementation
  in `src/materials/phasefield_smalldeform_nonAD/`).
- Parameter naming: snake_case, same style as `perm_exponent`, `wc`.
- Input-file idiom: hyphen-free keywords, SI units, explicit defaults.

### Numerical constraints
- Aperture must be non-negative; use the Macaulay bracket on `ε_n` since negative
  normal strain means compression/closure (no enhanced permeability).
- Cubic-law permeability uses `w²/12` (this is the 2D channel-flow result used in the
  existing code — preserve for consistency).
- Intrinsic permeability `k₀` is the floor: the enhanced permeability must satisfy
  `k_eff ≥ k₀` in every direction.
- Spurious anisotropy is a risk: when `ε_n` comes from a nearly isotropic strain
  tensor, the dominant eigenvector is ill-defined. Gate the enhancement on a
  damage threshold so the aperture is used only where `d > d_threshold` (default
  `d_threshold = 0` to preserve current behavior, user-tunable).

## Governing Equations

The formulation follows Heider (2021, Eng. Fract. Mech. 253:107881), Section 2.5,
equations (46)–(48), which are the standard Poiseuille-type, damage-gradient-based
smeared-crack permeability used by Miehe & Mauthe (2016), Wilson & Landis,
Heider & Sun, and others cited in the review.

### 1. Crack normal from the damage gradient

The unit normal to the smeared fracture surface at each quadrature point is

```
n_d = ∇d^S / |∇d^S|         (Heider eq. 46)
```

where `d^S` is the phase-field / damage variable. At undamaged points, `|∇d|` is
numerically zero; guard against division by zero with a tolerance `ε_grad`:

```
if |∇d| < ε_grad:    n_d = 0,   skip enhancement (fall back to K = K_poro)
else:                n_d = ∇d / |∇d|
```

**Alternative** (fallback): use the most-tensile principal-strain eigenvector
from the existing `_crack_rotation[_qp]` column 0. This is algebraically what
the current code already computes, but is less stable than `∇d / |∇d|` in
regions where the strain tensor is near-isotropic. Expose both via a parameter
`crack_normal_source = damage_gradient | principal_strain`, default
`damage_gradient`.

### 2. Normal strain and aperture (Heider eq. 47)

The normal strain to the fracture plane is

```
ε_nn = n_d · ε · n_d
```

Heider eq. (47) gives the opened-crack width as

```
w_c = ‖ h_c · (1 + n_d · ε · n_d) ‖ = h_c · |1 + ε_nn|       (Heider eq. 47)
```

where `h_c` is the initial length of a 1-D line element normal to the crack path.
For small tensile strains, `1 + ε_nn ≈ 1`, so the literal paper formula gives
`w_c ≈ h_c`, with damage localization provided by the `(d^S)^b` weighting in
eq. (48). This is the mesh-dependent "smeared-element aperture" interpretation.

A Macaulay-only variant (`w_c = h_c · <ε_nn>`, corresponding to displacement-jump
formulations Heider eqs. 53–55) was considered but rejected: it double-penalizes
intact material because `(d^S)^b` already suppresses fracture permeability at low
damage, so multiplying by `<ε_nn>` as well produces apertures ~10⁶× smaller
than the paper intends and makes the Heaviside gate `χ_d` redundant. Only the
Heider eq. (47) form is implemented in this plan.

The `h_c` source is user-selectable (matching the paper's "line element" notion):
- `characteristic_length_type = element_size` → `h_c = h` from
  `ElementLengthAux`-computed `mesh_size` aux variable (already present in
  `elasticity_E1d25.i:191–194, 328–333`). **This is the paper's intent** (the
  "1-D line element with initial length h_c").
- `characteristic_length_type = regularization_length` → `h_c = l` (MOOSE
  material property from `fracture_E1d25.i`; already propagated via `cli_args`).
- `characteristic_length_type = constant` → `h_c = characteristic_length_value`.

Default: `element_size` to match the paper.

### 3. Roughness-corrected fracture aperture (Heider eq. 46)

```
w_h = f_c · w_c · χ_d          (open crack, the case we model)
χ_d = H(d^S − d_threshold)     Heaviside step (1 if d ≥ threshold, else 0)
```

where `f_c ∈ (0, 1]` is a roughness correction factor (default `1.0`;
`correction_factor_fc` in the input) and `d_threshold` defaults to `0.5`
per Heider eq. (46). Expose `d_threshold` via parameter
`damage_threshold_for_permeability` (already in Phase 1 requirement 4). The
paper's closed-crack branch (`w_r`, residual aperture) is not needed for the
monotonic-damage pulse-power simulation; omit and document in the source.

### 4. Fracture permeability tensor (Heider eq. 46)

```
K_frac = (w_h² / 12) · (I − n_d ⊗ n_d)          (Heider eq. 46)
```

This is the cubic-law channel-flow permeability projected onto the **tangential
plane** of the fracture. `n_d ⊗ n_d` is the outer product rank-2 tensor of the
unit normal, and `I − n_d ⊗ n_d` is the tangential projector. Flow normal to
the fracture is not enhanced; flow parallel to the fracture is enhanced by the
cubic-law factor `w_h² / 12`. This replaces the rotation-matrix form I used in
the previous draft (the two are mathematically equivalent but the projector form
matches the paper and avoids building a full rotation tensor).

### 5. Total permeability (Heider eq. 48)

```
K = K_poro + (d^S)^b · K_frac                     (Heider eq. 48)
```

with `b ∈ {1, 2}` (paper: `b = 1` linear transition, `b = 2` quadratic
transition). Use `b = perm_exponent` in the input (matches the existing
parameter name), default kept at its current value `10` from
`elasticity_E1d25.i:42`. **Note**: the paper uses small integer `b`; the
existing code uses `b = 10`. The Phase 2 input-file update should decide
whether to stay at 10 (aggressive localization) or reduce to 2 (paper default).
Flag this as a user decision in Phase 2.

`K_poro` is the intrinsic matrix permeability `k₀ · I` (default `5e-19 I`),
optionally deformation-dependent (not implemented here; mentioned in paper
eq. 48 via `K_poro(d^S, n^S)` but using a constant matrix permeability is the
common simplification).

### 6. Simplified isotropic fallback

To allow regression runs that isolate the effect of anisotropy from the effect
of the new aperture, keep the isotropic fallback (identity tensor instead of
tangential projector):

```
K = K_poro + (d^S)^b · (w_h²/12) · I           (permeability_anisotropic = false)
```

### 7. Summary of changes relative to prior draft

| Item | Prior draft | Updated (matches Heider) |
| ---- | --- | --- |
| Crack normal | principal-strain eigvec | `∇d / |∇d|` (primary), eigvec (fallback) |
| Aperture | `<ε_n> · l_c` | `h_c · |1 + ε_nn|` (Heider eq. 47) |
| Anisotropy | rotation matrix `R · diag · Rᵀ` | tangential projector `I − n ⊗ n` |
| Blending | `k = k₀ + d^b(k_w − k₀)` | `K = K_poro + d^b · K_frac` (additive, paper eq. 48) |
| Roughness | not included | `f_c ∈ (0,1]` (default 1) |
| Damage gate | `d > threshold` (default 0) | `χ_d = H(d − 0.5)` (default 0.5) |
| Default `h_c` | regularization length | element size (matches paper's "line element") |

## Phase 1: Core Material — Normal-Strain Permeability Option

### Goal
`NDSmallDeformationIsotropicElasticity` exposes a third, physically grounded
permeability model that uses `<ε_n> · l_c` as the smeared-crack aperture and
optionally projects the enhanced permeability onto the two in-plane (tangential)
fracture directions, producing an anisotropic `effective_perm` tensor.

### Files to Create
- *(none)* — the change is contained in the existing material.

### Files to Modify
- `/Users/chunhuizhao/projects/farms_cdms/include/materials/phasefield_smalldeform_nonAD/NDSmallDeformationIsotropicElasticity.h`
  — add new params (bool + enum + real + real), change signature of
  `updatePermeabilityForCracking` to accept the principal-strain vector,
  add new member fields.
- `/Users/chunhuizhao/projects/farms_cdms/src/materials/phasefield_smalldeform_nonAD/NDSmallDeformationIsotropicElasticity.C`
  — register new `validParams` entries, initialize new members, validate
  mutual exclusivity, implement the new code branch inside
  `updatePermeabilityForCracking`.
- `/Users/chunhuizhao/projects/farms_cdms/include/materials/nonlocaldamage/FarmsComputeSmearedCrackingStressGradsSpectralSmallStrain.h`
  + corresponding `.C` in `src/materials/nonlocaldamage/` —
  **only if the input files in this run touch that material (grep first)**.
  If not used by the current cases, this file is *out of scope for Phase 1*; note
  in the risk section and revisit in a follow-up if a nonlocal variant is
  needed.

### Detailed Requirements

1. **Add a `PermeabilityModel` enum** in the header, inside the private section
   of `NDSmallDeformationIsotropicElasticity`:
   ```cpp
   enum class PermeabilityModel { none, exponential, darcy_poiseuille, normal_strain };
   const PermeabilityModel _permeability_model;
   ```
   Replace the three individual bool flags (`_exponential_permeability_model`,
   `_darcy_poiseuille_permeability_model`, and the new one) with this enum
   as the canonical storage. Keep the existing input-file bool parameters
   intact for backward compatibility; map them to the enum in the constructor.

2. **Add `MooseEnum` parameter** (`permeability_model`) in `validParams()`:
   ```cpp
   params.addParam<MooseEnum>(
       "permeability_model",
       MooseEnum("none exponential darcy_poiseuille normal_strain", "none"),
       "Permeability enhancement model to apply when porous_flow_coupling=true");
   ```
   If this new parameter is *not* `"none"`, it overrides the legacy bools with a
   `paramError` if they are inconsistent (bools set to `true` but enum is
   `"none"`, or vice versa).

3. **Add `crack_normal_source` MooseEnum parameter** (primary = damage gradient,
   per Heider eq. 46):
   ```cpp
   params.addParam<MooseEnum>(
       "crack_normal_source",
       MooseEnum("damage_gradient principal_strain", "damage_gradient"),
       "Source of the unit crack normal n_d. 'damage_gradient' uses grad(d)/|grad(d)| "
       "(Heider 2021 eq. 46). 'principal_strain' uses the eigenvector of the "
       "most-tensile principal strain (existing code path).");
   params.addParam<Real>(
       "damage_gradient_tolerance", 1e-30,
       "When |grad(d)| is below this tolerance, fall back to K = K_poro "
       "(no enhancement). Avoids division by zero in undamaged regions.");
   ```
   When `damage_gradient` is selected, read the gradient of the coupled
   phase-field variable via `coupledGradient("phase_field")`. Store as a
   new member `const VariableGradient & _grad_d;` and use
   `_grad_d[_qp]` inside `updatePermeabilityForCracking`.

4. **Add `characteristic_length_type` MooseEnum parameter**:
   ```cpp
   params.addParam<MooseEnum>(
       "characteristic_length_type",
       MooseEnum("element_size regularization_length constant", "element_size"),
       "Source of h_c in w_c = h_c*(...). 'element_size' matches the paper's "
       "1-D line element (default). 'regularization_length' uses l. 'constant' "
       "uses a user-supplied value.");
   params.addParam<MaterialPropertyName>(
       "regularization_length_name", "l",
       "Name of the material property holding the regularization length "
       "(used when characteristic_length_type = regularization_length).");
   params.addCoupledVar(
       "element_size_variable",
       "Name of an aux variable holding the local element size h (used when "
       "characteristic_length_type = element_size).");
   params.addParam<Real>(
       "characteristic_length_value", -1.0,
       "Constant h_c (used when characteristic_length_type = constant).");
   ```

5. **Add anisotropy, threshold, roughness, and blending-exponent parameters**
   (Heider eqs. 46, 48):
   ```cpp
   params.addParam<bool>(
       "permeability_anisotropic", true,
       "If true, K_frac = (w^2/12)*(I - n_d (x) n_d) (Heider eq. 46). "
       "If false, K_frac = (w^2/12)*I (isotropic fallback).");
   params.addParam<Real>(
       "damage_threshold_for_permeability", 0.5,
       "Heaviside gate chi_d = H(d - threshold) (Heider eq. 46). "
       "K_frac is zero where d < threshold. Default 0.5 matches the paper.");
   params.addParam<Real>(
       "correction_factor_fc", 1.0,
       "Roughness correction factor f_c in w_h = f_c*w_c*chi_d (Heider eq. 46). "
       "Default 1.0 (smooth walls).");
   ```
   The damage-localization exponent `b` is the existing `perm_exponent` parameter;
   no new input is needed for it.

6. **Constructor initialization** (declare-and-initialize in the C++ init list):
   - Read `permeability_model` MooseEnum → `_permeability_model`.
   - Read `crack_normal_source` → `_normal_source`.
   - Read `characteristic_length_type` → `_lc_type`.
   - Bind `_grad_d = coupledGradient("phase_field")` (always, since the
     `phase_field` variable is already coupled; this is a zero-cost addition
     when `damage_gradient` is not selected).
   - Conditionally retrieve the regularization-length material property
     (`_l_mat_prop`), the coupled element-size variable (`_h_elem`), or the
     constant value (`_lc_const`).
   - Read `permeability_anisotropic` → `_perm_anisotropic`.
   - Read `damage_threshold_for_permeability` → `_d_perm_threshold`.
   - Read `correction_factor_fc` → `_fc`.
   - Read `damage_gradient_tolerance` → `_grad_d_tol`.

7. **Mutual exclusivity and validity checks** (in the constructor body):
   - If `_permeability_model == PermeabilityModel::normal_strain`:
     - Require `_porous_flow_coupling == true`, else `paramError`.
     - Require `_perm_exponent > 0`, else `paramError`.
     - If `_lc_type == element_size`: require that
       `isCoupled("element_size_variable")` is true, else `paramError`.
     - If `_lc_type == constant`: require `_lc_const > 0`, else `paramError`.
     - Require `_fc > 0 && _fc <= 1`, else `paramError`.
     - Require `_d_perm_threshold >= 0 && _d_perm_threshold < 1`, else
       `paramError`.
   - If `_permeability_model == PermeabilityModel::none` AND `_porous_flow_coupling`
     is true: preserve the existing `paramError` message ("no permeability model
     is selected"), extended to mention `permeability_model` as a valid
     alternative.

8. **Change signature of `updatePermeabilityForCracking`** to accept the
   principal strain vector (still needed for the `principal_strain` fallback
   path and for the existing Darcy–Poiseuille branch to remain unchanged):
   ```cpp
   // OLD: void updatePermeabilityForCracking();
   // NEW:
   void updatePermeabilityForCracking(const RealVectorValue & strain_in_crack_dir);
   ```
   Update the call site at `NDSmallDeformationIsotropicElasticity.C:250` to
   pass the local `strain_in_crack_dir` that was computed two lines above.

9. **Implement the `normal_strain` branch** (Heider eqs. 46–48) inside
    `updatePermeabilityForCracking`:
    ```cpp
    else if (_permeability_model == PermeabilityModel::normal_strain)
    {
      // (a) Determine unit crack normal n_d
      RealVectorValue n_d;
      bool have_normal = false;

      if (_normal_source == CrackNormalSource::damage_gradient)
      {
        const Real gnorm = _grad_d[_qp].norm();
        if (gnorm > _grad_d_tol)
        {
          n_d = _grad_d[_qp] / gnorm;
          have_normal = true;
        }
      }
      else   // principal_strain fallback
      {
        // _crack_rotation[_qp] column 0 is the most-tensile eigenvector
        // (see computeCrackStrainAndOrientation).
        n_d = _crack_rotation[_qp].column(0);
        have_normal = (n_d.norm() > 0.0);
      }

      const Real k0 = _intrinsic_permeability;
      const RankTwoTensor I = RankTwoTensor::Identity();
      const Real d = _d[_qp];

      // Heaviside gate chi_d (Heider eq. 46)
      const Real chi_d = (d >= _d_perm_threshold) ? 1.0 : 0.0;

      // If no valid normal OR below threshold OR d = 0: fall back to matrix perm
      if (!have_normal || chi_d == 0.0 || d <= 0.0)
      {
        _effective_perm[_qp] = k0 * I;   // K_poro only
        return;
      }

      // (b) Normal strain eps_nn = n_d . eps . n_d
      const Real eps_nn = _elastic_strain[_qp].contractionWithVector(n_d) * n_d;
      // NOTE: if `contractionWithVector` is not available, use the explicit
      // double-contraction:  eps_nn = n_d_i * eps_ij * n_d_j via two loops.

      // (c) Aperture w_c (Heider eq. 47): w_c = h_c * |1 + eps_nn|
      const Real h_c = getCharacteristicLength();
      const Real w_c = h_c * std::abs(1.0 + eps_nn);

      // (d) Roughness-corrected aperture (Heider eq. 46): w_h = f_c * w_c * chi_d
      const Real w_h = _fc * w_c * chi_d;
      const Real k_w = w_h * w_h / 12.0;

      // (e) Fracture permeability tensor
      RankTwoTensor K_frac;
      if (_perm_anisotropic)
      {
        // Tangential projector (Heider eq. 46): K_frac = k_w * (I - n_d (x) n_d)
        RankTwoTensor n_outer_n;
        for (unsigned int i = 0; i < 3; ++i)
          for (unsigned int j = 0; j < 3; ++j)
            n_outer_n(i, j) = n_d(i) * n_d(j);
        K_frac = k_w * (I - n_outer_n);
      }
      else
      {
        K_frac = k_w * I;   // isotropic fallback
      }

      // (f) Damage-weighted total perm (Heider eq. 48)
      //     K = K_poro + (d^b) * K_frac,  b = _perm_exponent
      const Real weight = std::pow(d, _perm_exponent);
      _effective_perm[_qp] = k0 * I + weight * K_frac;
    }
    ```
    (Adapt the `RealTensorValue` ↔ `RankTwoTensor` conversion used in the
    existing branches; the current code stores into `_effective_perm` from a
    `RankTwoTensor` via implicit conversion, keep that convention.)

10. **Add a private helper**:
    ```cpp
    Real getCharacteristicLength() const;   // .h
    Real NDSmallDeformationIsotropicElasticity::getCharacteristicLength() const
    {
      switch (_lc_type)
      {
        case LcType::regularization_length:
          return (*_l_mat_prop)[_qp];
        case LcType::element_size:
          return _h_elem[_qp];                 // coupledValue
        case LcType::constant:
          return _lc_const;
      }
      mooseError("Unknown characteristic_length_type.");
    }
    ```

11. **Do not modify the existing `exponential_permeability_model` or
    `darcy_poiseuille_permeability_model` branches**. They remain intact for
    regression comparability. The legacy bool parameters still work; when both a
    legacy bool and the new `permeability_model` enum are set, the enum wins and
    a warning is printed (`mooseWarning`).

### Interfaces

- Public (input file) interface adds these parameters:
  - `permeability_model = normal_strain | exponential | darcy_poiseuille | none`
  - `crack_normal_source = damage_gradient | principal_strain` (default `damage_gradient`)
  - `damage_gradient_tolerance` (Real, default `1e-30`)
  - `characteristic_length_type = element_size | regularization_length | constant`
    (default `element_size`)
  - `regularization_length_name` (MaterialPropertyName, default `"l"`)
  - `element_size_variable` (coupled AuxVariable name)
  - `characteristic_length_value` (Real)
  - `permeability_anisotropic` (bool, default `true`)
  - `damage_threshold_for_permeability` (Real, default `0.5`)
  - `correction_factor_fc` (Real, default `1.0`)
- Private C++ interface:
  - Scoped enums: `PermeabilityModel`, `CrackNormalSource`, `LcType`.
  - `updatePermeabilityForCracking(const RealVectorValue & strain_in_crack_dir)`
    — new signature (argument retained for `principal_strain` path and existing
    Darcy–Poiseuille branch).
  - `getCharacteristicLength()` — new private method.
  - New member: `const VariableGradient & _grad_d;`
- Material property `effective_perm` name, type, and semantics unchanged.

### Edge Cases to Handle
1. **Undamaged region (`|∇d| < tol`)**: damage gradient ill-defined. Fall
   back to `K = k₀ · I` (matrix permeability only). This naturally yields
   `K_frac = 0` and is consistent with `d = 0` case.
2. **Compressive normal strain (`ε_nn < 0`)**:
   the paper formula gives `w_c = h_c · |1 + ε_nn|`, which is
   well-defined for any `ε_nn > −1`. For physically realistic strain magnitudes
   (`|ε_nn| ≪ 1`), `w_c ≈ h_c`. No special handling needed beyond the absolute
   value in the formula.
3. **Strongly compressive strain (`ε_nn < −1`)** — unphysical in small-strain:
   the Heider formula gives `w_c = h_c · |1 + ε_nn| = h_c · |ε_nn − (−1)|`.
   If the simulation ever reaches this regime, the permeability would grow
   again due to the absolute value. Add a `mooseWarning` when `ε_nn < −0.5`
   to alert the user that small-strain assumptions are being violated; clamp
   `1 + ε_nn` to a floor of `0` (so `w_c = 0` below that).
4. **Zero damage (`d = 0`)**: `(d)^b = 0` → `K_frac` weight is zero → `K = k₀·I`. ✓
5. **`d` between 0 and threshold (e.g. 0.3 with threshold 0.5)**: `χ_d = 0` →
   `K_frac = 0` → `K = k₀·I`. ✓
6. **Fully damaged + near-zero normal strain**: `w_c ≈ h_c` (not zero), so
   `K_frac ≈ (h_c²/12) · (I − n_d ⊗ n_d)` is non-zero. This matches the
   paper's intent: a fully damaged element represents an open fracture of
   width `h_c`, with cubic-law permeability of that width.
7. **2D mesh with `LIBMESH_DIM == 3`**: in 2D the damage gradient has
   `∂d/∂z = 0`, so `n_d` lies in the xy-plane. The tangential projector
   `I − n_d ⊗ n_d` still has a contribution along `e_z` (the out-of-plane
   direction), which correctly represents flow parallel to the crack in 3D.
   **Verify by test** — the in-plane entries `K_xx`, `K_yy`, `K_xy` should
   match the 2D analytical calculation.
8. **Negative `_perm_exponent`** (legacy default `-1.0` if not set): already
   validated to be `>0` by the existing parameter checks when the
   Darcy–Poiseuille model is used; replicate that validation when the
   normal-strain model is selected.
9. **Input file omits `porous_flow_coupling = true` but sets
   `permeability_model = normal_strain`**: error out with a clear message.
10. **`characteristic_length_type = element_size` but `mesh_size` aux var is
    0**: first timestep may not have computed `mesh_size` yet. Guard the
    computation by ordering `ElementLengthAux` at `TIMESTEP_BEGIN` (already
    done in `elasticity_E1d25.i:332`) and at `INITIAL` (add this if not
    already present).

### Acceptance Criteria
- [ ] `NDSmallDeformationIsotropicElasticity` compiles after the header/source
      changes.
- [ ] Existing `exponential_permeability_model` and `darcy_poiseuille_permeability_model`
      test runs (if any exist in `test/tests/.../`) produce byte-identical
      `effective_perm` values before and after the refactor. Specifically:
      grep `permeability` in `/Users/chunhuizhao/projects/farms_cdms/test/`
      and re-run every match through `run_tests`.
- [ ] A minimal unit-style regression input using `permeability_model = normal_strain`,
      `crack_normal_source = damage_gradient`, `permeability_anisotropic = true`,
      a hand-set linear damage profile `d(x,y) = x` (so `∇d = e_x` and
      `n_d = e_x`), uniaxial strain `ε_xx = ε_set`, and damage `d = 0.5` at
      the QP yields:
      `K_xx ≈ k₀` (normal direction killed by tangential projector),
      `K_yy = K_zz ≈ k₀ + 0.5^b · (f_c · h_c · |1 + ε_set|)² / 12`,
      `K_ij = 0` for `i ≠ j`.
- [ ] Same unit test with `permeability_anisotropic = false` yields
      `K_xx = K_yy = K_zz` (isotropic fallback).
- [ ] A rotation test with `∇d` at 45° in the xy-plane yields the correct
      off-diagonal `K_xy`. Specifically, with `n_d = (1,1,0)/√2` and
      `k_frac_scalar ≡ (d^b)·k_w`, the expected tensor is
      `K = k₀·I + k_frac_scalar · (I − n_d ⊗ n_d)`, which evaluates to
      `K_xx = K_yy = k₀ + k_frac_scalar/2`, `K_xy = −k_frac_scalar/2`,
      `K_zz = k₀ + k_frac_scalar`.
- [ ] Running `elasticity_E1d25.i` + `fracture_E1d25.i` with the input file
      set to the legacy `darcy_poiseuille_permeability_model` (i.e. without
      the Phase 2 input changes) produces a regression that matches the
      current CSV `full_energy` within numerical noise (≤ 1e-10 rel. diff),
      proving that the legacy code path is unaffected by the refactor.

### Dependencies
- Depends on: *(nothing — this is a self-contained material extension).*
- Required by: Phase 2, Phase 3.

## Phase 2: Wire the New Model into the Pulse-Power Input Files

### Goal
`elasticity_E1d25.i` and `fracture_E1d25.i` drive the coupled simulation using
the new normal-strain aperture model, producing a physically grounded
permeability field that does not require the phenomenological `wc` value.

### Files to Modify
- `/Users/chunhuizhao/projects/farms_cdms/pulsepower/pf_code2d_porousflow/parametric_study/case_cf1_domain1x/elasticity_E1d25.i`
- `/Users/chunhuizhao/projects/farms_cdms/pulsepower/pf_code2d_porousflow/parametric_study/case_cf1_domain1x/fracture_E1d25.i`
  — **inspect but do not modify**; `l` is already propagated via
  `cli_args = 'Gc_const=${Gc_const};l=${l}'` (line 56 of `elasticity_E1d25.i`),
  so no change needed here unless tests demand one.

### Detailed Requirements

1. **Replace** the block in `elasticity_E1d25.i` at lines 660–666:
   ```ini
   ##-----darcy_poiseuille_permeability_model-----##
   darcy_poiseuille_permeability_model = true
   intrinsic_permeability = ${intrinsic_permeability}
   wc = ${wc}
   perm_exponent = ${perm_exponent}
   ```
   with the Heider-2021 normal-strain formulation:
   ```ini
   ##-----normal_strain_permeability_model (Heider 2021, eqs. 46-48)-----##
   permeability_model = normal_strain
   intrinsic_permeability = ${intrinsic_permeability}
   perm_exponent = ${perm_exponent}          # exponent b in K = K_poro + d^b*K_frac
   crack_normal_source = damage_gradient     # n_d = grad(d)/|grad(d)|
   characteristic_length_type = element_size # h_c = element size (paper default)
   element_size_variable = mesh_size         # reuse existing mesh_size AuxVariable
   permeability_anisotropic = true           # K_frac = (w^2/12)(I - n_d (x) n_d)
   damage_threshold_for_permeability = 0.5   # chi_d = H(d - 0.5) per eq. (46)
   correction_factor_fc = 1.0                # smooth-walled default
   ```

2. **User-facing decision on `perm_exponent`**: the existing value is `10`
   (aggressive localization). Heider eq. (48) uses `b ∈ {1, 2}` (linear /
   quadratic). Document both options in the file with a comment:
   ```ini
   # perm_exponent = 10 # aggressive localization (existing pulse-power value)
   # perm_exponent = 2  # quadratic (Heider 2021 eq. 48)
   # perm_exponent = 1  # linear    (Heider 2021 eq. 48)
   perm_exponent = 10
   ```
   Keep the default at `10` to preserve comparability with existing runs;
   note that switching to `1` or `2` is a separate calibration study.

3. **Remove now-unused top-of-file variables** `wc` and the comment referencing
   the Darcy–Poiseuille ultimate crack opening (lines 40–42 of
   `elasticity_E1d25.i`). Replace with a comment documenting the new model:
   ```ini
   ##Heider-2021 normal-strain permeability model:
   ##  aperture w_c = h_c * |1 + n_d . eps . n_d|,   K_frac = (w_c^2/12) (I - n_d (x) n_d)
   ##  total K = k0*I + d^b * K_frac
   perm_exponent = 10 # damage localization exponent b (see below)
   ```

4. **Do not** remove `intrinsic_permeability` — it is `k₀` in `K_poro = k₀·I`.

5. **Verify `mesh_size` AuxVariable availability at `INITIAL`**: the current
   input defines `[./max]` (`ElementLengthAux`) with
   `execute_on = TIMESTEP_BEGIN` only (line 332). Add `INITIAL` to cover the
   first material evaluation:
   ```ini
   execute_on = 'INITIAL TIMESTEP_BEGIN'
   ```

6. **`fracture_E1d25.i` requires no change.** Verify via `diff` that only the
   `elasticity_E1d25.i` file is modified in this phase.

7. **Add a short header comment** at the top of the modified file (lines 1–10)
   documenting the change:
   ```ini
   # Permeability enhancement: Heider (2021) normal-strain formulation (eqs. 46-48)
   # Crack normal n_d = grad(d)/|grad(d)|, aperture w_c = h_c*|1 + n_d.eps.n_d|,
   # fracture perm K_frac = (w_c^2/12)(I - n_d (x) n_d), total K = k0*I + d^b*K_frac.
   # Replaces the prior simplified w = d*wc formulation.
   ```

### Interfaces
- No new MOOSE-level interfaces; this is a pure input-file update that consumes
  the Phase 1 interfaces.

### Edge Cases to Handle
- If the user runs `elasticity_E1d25.i` with an old `farms` binary that does not
  include Phase 1, the parse will fail with an "unknown parameter
  `permeability_model`" error — accept this as the required failure mode.
- The static-solve output file `static_solve_out.e` (loaded in
  `[UserObjects]/init_sol_components`) was generated with the old
  permeability model. This is fine: the static solve only supplies initial
  conditions for `disp_x`, `disp_y`, `pp`, and elastic-strain components.
  **However**, if the static solve itself also uses `elasticity_E1d25.i`-style
  permeability parameters, the static input (`static_solve.i`) must be updated
  consistently before regeneration. Grep `static_solve.i` for
  `darcy_poiseuille_permeability_model` before finishing this phase and update
  if present.

### Acceptance Criteria
- [ ] `moose-opt -i elasticity_E1d25.i --mesh-only` parses without error.
- [ ] Running the first 10 time steps (`end_time = 1e-7`, `dt = 1e-8`) to
      convergence produces finite values in `effective_perm00_aux`,
      `effective_perm11_aux`, `effective_perm01_aux` for every quadrature point
      (no NaNs or infs).
- [ ] Qualitative comparison against the prior (`w = d · wc`) run:
      (a) where `d < 0.5` (below the Heaviside gate), the new permeability
      should be exactly `k₀` everywhere — stricter localization than the old
      model, which applied `d^n · (...)` even at small `d`.
      (b) where `d ≥ 0.5`, the new in-plane permeability should be anisotropic
      (unequal `effective_perm00_aux` vs `effective_perm11_aux` at points
      where the damage gradient is not aligned with x or y), whereas the old
      model was isotropic. This is the signature test that the tangential
      projector is active.

### Dependencies
- Depends on: Phase 1.
- Required by: Phase 3.

## Phase 3: Unit Tests and Regression

### Goal
Lock the new code path behind automated tests so that future refactors do not
silently break the normal-strain formulation.

### Files to Create
- `/Users/chunhuizhao/projects/farms_cdms/test/tests/materials/permeability_normal_strain/isotropic_axis.i`
- `/Users/chunhuizhao/projects/farms_cdms/test/tests/materials/permeability_normal_strain/anisotropic_axis.i`
- `/Users/chunhuizhao/projects/farms_cdms/test/tests/materials/permeability_normal_strain/anisotropic_rotated.i`
  (the critical 45°-rotation test)
- `/Users/chunhuizhao/projects/farms_cdms/test/tests/materials/permeability_normal_strain/principal_strain_fallback.i`
- `/Users/chunhuizhao/projects/farms_cdms/test/tests/materials/permeability_normal_strain/legacy_darcy_poiseuille.i`
- `/Users/chunhuizhao/projects/farms_cdms/test/tests/materials/permeability_normal_strain/tests`
  (MOOSE `tests` spec file listing all five tests with their gold files)
- Corresponding `gold/*.e` files (generated with the Phase 1 implementation,
  hand-verified against the analytical formulas in the requirements, and
  checked into the repository).

### Files to Modify
- *(none)* — tests are self-contained.

### Detailed Requirements

All unit tests use a single QUAD4 element of size 1 × 1 m and prescribe the
damage field `d` via a `FunctionAux` so the damage gradient `∇d` is
analytically known at the QP. Strain is prescribed via `FunctionDirichletBC`
on all four sides.

Common material parameters (unless overridden per test):
  - `intrinsic_permeability = 5e-19`
  - `perm_exponent = 2`             (quadratic, per Heider eq. 48)
  - `damage_threshold_for_permeability = 0.5`
  - `correction_factor_fc = 1.0`
  - `crack_normal_source = damage_gradient`
  - `characteristic_length_type = constant`
  - `characteristic_length_value = 1e-3` (so `h_c = 1e-3` m)
  - `permeability_model = normal_strain`

1. **Axis-aligned damage gradient, isotropic fallback (`isotropic_axis.i`)**:
   - Damage profile: `d(x, y) = 0.5 + 0.4·x` (on [0,1] element → `d` from 0.5
     to 0.9). `∇d = (0.4, 0, 0)` → `n_d = (1, 0, 0)`.
   - Strain: `disp_x = ε_set · x`, `disp_y = 0`, `ε_set = 1e-3`.
   - `ε_nn = n_d · ε · n_d = ε_xx = 1e-3`.
   - `permeability_anisotropic = false`.
   - Expected at any QP (using QP damage value `d_qp`):
     - `w_c = h_c · |1 + ε_nn| = 1e-3 · 1.001 = 1.001e-3` m.
     - `w_h = f_c · w_c · χ_d = w_c · 1 = 1.001e-3` (since `d_qp ≥ 0.5`).
     - `k_w = w_h²/12 ≈ 8.35e-8` m².
     - `K_frac = k_w · I`.
     - `K = k₀·I + d_qp² · k_w · I = (5e-19 + d_qp² · 8.35e-8) · I`.
   - Assertion: `MaterialRealTensorValueAux` entries match the analytical
     formula at the QP within absolute tolerance `1e-12` m² (note: the
     magnitudes here are large because `h_c = 1e-3` — this is intentional
     so exodiff tolerances don't need to be near machine epsilon).

2. **Axis-aligned damage gradient, anisotropic, uniaxial strain
   (`anisotropic_axis.i`)**:
   - Same damage and strain as Test 1.
   - `permeability_anisotropic = true`.
   - `n_d = (1, 0, 0)`, so the tangential projector is
     `I − n_d ⊗ n_d = diag(0, 1, 1)`.
   - Expected at QP:
     - `K_xx = k₀` (normal to crack, no enhancement).
     - `K_yy = K_zz = k₀ + d_qp² · k_w`.
     - `K_xy = K_xz = K_yz = 0`.
   - Assertion: exodiff against analytical gold.

3. **Rotation test — 45° damage gradient (`anisotropic_rotated.i`)**:
   - Damage profile: `d(x, y) = 0.5 + 0.2·(x + y)` →
     `∇d = (0.2, 0.2, 0)` → `n_d = (1, 1, 0)/√2`.
   - Strain: `disp_x = 0.5·ε_set·(x + y)`, `disp_y = 0.5·ε_set·(x + y)` →
     `ε_xx = ε_yy = ε_xy = 0.5·ε_set`, with `ε_set = 2e-3`.
   - `ε_nn = n_d · ε · n_d = 0.5·ε_set·(1+1+2·1)/2 = ε_set = 2e-3`.
     (Compute explicitly: `n_d · ε · n_d = (1/2)(ε_xx + 2·ε_xy + ε_yy) = ε_set`.)
   - `w_c = h_c · |1 + ε_nn| = 1e-3 · 1.002 = 1.002e-3`.
   - `k_w = w_c²/12 ≈ 8.367e-8`.
   - Tangential projector with `n_d = (1,1,0)/√2`:
     - `(n_d ⊗ n_d)` = `[[1/2, 1/2, 0], [1/2, 1/2, 0], [0, 0, 0]]`.
     - `I − n_d ⊗ n_d` = `[[1/2, −1/2, 0], [−1/2, 1/2, 0], [0, 0, 1]]`.
   - Let `α = d_qp² · k_w`. Expected:
     - `K_xx = K_yy = k₀ + α/2`.
     - `K_xy = −α/2`.
     - `K_zz = k₀ + α`.
     - `K_xz = K_yz = 0`.
   - **This test is the critical one for verifying the tangential projector
     is implemented correctly in the global frame.**
   - Assertion: exodiff against analytical gold, tolerance `1e-14 · max(|K|)`
     relative.

4. **Principal-strain fallback (`principal_strain_fallback.i`)**:
   - Uniform damage `d = 0.7` (constant → `∇d = 0`).
   - Strain: `ε_xx = 1e-3`, all others zero.
   - `crack_normal_source = principal_strain`.
   - `|∇d| = 0` < tolerance → under `damage_gradient` source this would fall
     back to matrix perm. But with `principal_strain` source, the max-tensile
     eigenvector is `e_x`, so `n_d = (1, 0, 0)` and the test mirrors Test 2.
   - Assertion: same expected tensor as Test 2 (with `d_qp = 0.7` constant).

5. **Legacy regression (`legacy_darcy_poiseuille.i`)**:
   - A single-element test that exercises the old `darcy_poiseuille_permeability_model`
     branch with `wc = 1e-6`, `d = 0.5`, `perm_exponent = 10`. Expected
     `effective_perm` from the original formula
     `k_eff = k₀ + d^n · ((d·wc)²/12 − k₀)`, then rotated — but since the
     rotation is of an isotropic tensor it's a no-op.
   - Expected: `k_eff = (5e-19 + 0.5^10 · ((0.5·1e-6)²/12 − 5e-19)) · I`.
   - This test must pass **both before and after the refactor** to prove the
     legacy code path is preserved.

### Acceptance Criteria
- [ ] `./run_tests --re permeability_normal_strain` passes.
- [ ] The legacy Darcy–Poiseuille test (either preexisting or newly added) passes
      with identical gold output to its pre-refactor version.
- [ ] All three hand-computed analytical permeability values (isotropic,
      anisotropic, 45°-rotated) agree with the simulation output within the
      stated tolerance.

### Dependencies
- Depends on: Phase 1 (code changes), Phase 2 optional (not required).
- Required by: *(nothing — this is the final phase).*

## Testing Strategy

1. **Phase 1 code review**:
   - Compile with `-Wall -Werror`.
   - Run `./run_tests -j8` across `test/tests/materials/` to confirm no existing
     material tests regress.
   - Hand-verify that `effective_perm` before and after the refactor is
     byte-identical for the legacy Darcy–Poiseuille branch by running the
     existing `elasticity_E1d25.i` *with the old input block* — it should
     produce the same CSV output as the pre-refactor binary.

2. **Phase 2 smoke test**:
   - Run `elasticity_E1d25.i` for 10 steps (short circuit the `end_time` to
     `1e-7`) and dump `effective_perm00_aux`, `effective_perm11_aux`, and
     `effective_perm01_aux`. Checks:
     - Below the Heaviside gate (`d < 0.5`): all three should equal `k₀`
       exactly (no enhancement, matrix perm only).
     - Above the gate (`d ≥ 0.5`): the in-plane components should follow
       `k_ij = k₀·δ_ij + d^b · (w_h²/12) · (δ_ij − n_i n_j)`.
     - The off-diagonal `effective_perm01_aux` should be **non-zero** wherever
       the damage gradient is not aligned with x or y — this is the signature
       that the anisotropic tangential projector is active.
     - At quadrature points where `∇d = 0` (truly uniform damage region):
       the fallback sets `K = k₀·I` (off-diagonals zero). Check this in the
       far-field (away from the damage localization).

3. **Phase 3 unit tests**:
   - Exact analytical gold values computed by hand, compared via `exodiff`
     with absolute tolerance chosen per test (generally `1e-14` relative to
     the largest expected component).

4. **Energy-balance sanity** (Phase 2 full run):
   - The CSV postprocessors `full_energy`, `full_input_energy`, etc. should
     remain in approximate balance (drift < 1% over the full run). A
     catastrophic balance failure indicates a unit-error or sign-error in the
     permeability tensor (e.g. tangential projector built with `+ n⊗n` instead
     of `− n⊗n`, or the damage weight `d^b` applied to the matrix perm instead
     of the fracture perm).

## Risk Assessment

1. **Tangential-projector sign/ordering error**. With `n_d` constructed either
   from the damage gradient or from the principal-strain eigenvector, the
   tangential projector `I − n_d ⊗ n_d` must correctly zero the normal
   direction when the projector is applied to a scalar-weighted identity.
   An off-by-one indexing error in the outer product would break the rotation
   test. **Detection**: the 45°-rotation unit test (Test 4 in Phase 3) is
   designed specifically for this — the expected `K_xy = −α/2` has a unique
   sign and magnitude that any incorrect projector would violate.
   **Mitigation**: the 45°-rotation test is mandatory.

2. **Damage gradient is noisy at low damage**. Near the fracture front where
   `d` transitions from ~0 to ~1, `∇d` can be poorly conditioned (numerical
   spikes). **Detection**: visualize `n_d_x`, `n_d_y` aux components in the
   first run. **Mitigation**: the Heaviside gate `χ_d = H(d − 0.5)` already
   suppresses enhancement below `d = 0.5`, and the `damage_gradient_tolerance`
   parameter kills the enhancement when `|∇d| < tol`. If spurious anisotropy
   persists above `d = 0.5`, raise the threshold to `0.7` or switch to
   `crack_normal_source = principal_strain` (less noisy in some regimes).

3. **Paper formula gives non-zero aperture at zero strain**.
   The literal paper formula `w_c = h_c · |1 + ε_nn|` evaluates to `h_c`
   when `ε_nn = 0`, which means even a newly damaged element (where `d ≥ 0.5`
   but strain is tiny) will have `K_frac = (h_c)²/12 · (I − n_d ⊗ n_d)`.
   For pulse-power simulations with `h_c ~ 10⁻⁴ m`, this corresponds to
   `k_w ~ 8.3 × 10⁻¹⁰ m²`, about 9 orders of magnitude above the matrix
   permeability. **This is the paper's intended behavior** (a fully damaged
   element represents an open fracture whose width equals the element size),
   but is also orders of magnitude larger than what the current `w = d·wc`
   formulation produces. **Mitigation**: the damage weight `(d^S)^b`
   (`perm_exponent`) and the Heaviside gate `χ_d` together limit the
   enhancement to regions where `d ≥ 0.5`. If the resulting permeabilities
   destabilize the solver or produce nonphysical pressure evolution, the
   primary tuning knobs are `perm_exponent` (try `b = 1` or `b = 2` per
   Heider eq. 48) and `correction_factor_fc` (physically `<1` for rough
   crack walls). Document these tradeoffs in the input-file comment.

4. **Characteristic-length choice is physics-loaded**. The paper's "1-D line
   element with initial length h_c" strongly implies `h_c = element size`,
   making the formulation mesh-dependent (finer mesh → smaller aperture →
   lower fracture permeability). This can be desirable (it ties the smeared
   crack's permeability to the mesh-resolved crack band) or undesirable (it
   breaks mesh-convergence studies). **Mitigation**: default to `element_size`
   (paper's choice) but expose `regularization_length` and `constant` for
   users who want mesh-independent permeability. Document this explicitly.

5. **Static solve uses legacy permeability**. The initial conditions are loaded
   from `static_solve_out.e`. If the static solve input (`static_solve.i`) is
   not updated consistently, the transient run will inherit a stress/pressure
   state that was computed under a different permeability law. **Detection**:
   grep `static_solve.i` for `darcy_poiseuille` before Phase 2 completion.
   **Mitigation**: Phase 2 edge-case 2 explicitly calls this out.

6. **Nonlocal-damage variant (`FarmsComputeSmearedCrackingStressGradsSpectralSmallStrain`)
   has its own copy of the permeability code**
   (`src/materials/nonlocaldamage/FarmsComputeSmearedCrackingStressGradsSpectralSmallStrain.C:344–367`).
   If the user later switches the simulation to the nonlocal variant, the new
   normal-strain branch will not be available. **Mitigation**: out of scope for
   this plan; documented in the header comment of the modified source file so
   that the next engineer sees it. Revisit if the pulse-power cases adopt the
   nonlocal variant.

7. **Existing `_wc` and legacy `darcy_poiseuille_permeability_model` flag remain
   after the refactor**. They introduce dead-code branches if no input uses
   them, but removing them would break backward compatibility with the
   sibling input files in `parametric_study/case_cf1_domain1x/` (e.g.
   `elasticity_E1d25_pulse2em5.i` etc., which may still use Darcy–Poiseuille).
   **Mitigation**: leave them in place; do not touch the sibling inputs in
   Phase 2.

8. **`perm_exponent = 10` is not from Heider**. Heider eq. (48) shows
   `b ∈ {1, 2}`. The existing pulse-power value `10` produces extremely
   localized enhancement (at `d = 0.9`, weight is `0.9^10 ≈ 0.35`; at
   `d = 0.99`, weight is `0.99^10 ≈ 0.90`). This may be physically excessive.
   **Mitigation**: Phase 2 requirement 2 flags the exponent choice for user
   review. A calibration study (changing `b` to `1` or `2` and comparing
   field evolution) is listed as out-of-scope here but recommended as a
   follow-up.

## Out-of-Scope / Investigation Items

- **Reference verified**: Heider (2021) "A review on phase-field modeling of
  hydraulic fracturing", Engineering Fracture Mechanics 253:107881, Section 2.5
  (equations 46–48) is the authoritative formulation used in this plan.
  Key references cited therein for further detail:
  - **Heider & Sun (2020)** [ref 82 in the paper] — the 1-D line element
    implementation for `w_c`, which is what this plan most closely follows.
  - **Miehe, Mauthe (2016)** [ref 31] — original spectral-decomposition
    approach; fully saturated poroelasticity + PFM.
  - **Wilson & Landis** [refs 63, 65, 66] — Ginzburg–Landau-type formulation,
    same family of permeability enhancement.
  Reading the Heider & Sun (2020) paper directly is the recommended follow-up
  if fine-tuning of `f_c`, `χ_d`, or the closed-crack residual aperture `w_r`
  becomes relevant.

- **Calibration of `perm_exponent` (`b`)**: the Heider paper uses `b = 1` or
  `b = 2`. The existing pulse-power input uses `b = 10`. A sensitivity study
  on this parameter is recommended as a follow-up to the main implementation,
  independent of the code change itself.

- **Closed-crack residual aperture `w_r`** (Heider eq. 46, second branch)
  is omitted from this plan because the pulse-power simulation is in a
  tension-dominated, monotonically damaging regime. If the user extends to
  cyclic loading or crack closure, `w_r` and its companion Heaviside
  `(f_c · w_r) · χ_d` for closed cracks must be added.

- **Nonlocal-damage material** (`FarmsComputeSmearedCrackingStressGradsSpectralSmallStrain`)
  is left untouched; if the user's target workflow will switch to nonlocal
  damage, a follow-up plan must mirror Phase 1 in that material.
