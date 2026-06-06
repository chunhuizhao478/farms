# Permeability-formula comparison set: 5 cases × 2 meshes

This directory holds a controlled comparison of the `normal_strain` fracture-
permeability model in the 2D hydro-mechanical pulse problem. There are **five
case directories**, each run on **two meshes** (10 runs total). Every case is
identical except for the one parameter under study, so differences in the
output fields are attributable to that parameter alone.

`baseline/` is the shared control (original model, no new features). **Group A**
varies the *crack-normal definition* `n_F` off the baseline; **Group B** varies
the *porosity-update law* `φ` off the baseline.

| Case dir                | Group | Crack normal `n_F`                                   | Porosity law `φ`                         |
|-------------------------|-------|-----------------------------------------------------|------------------------------------------|
| `baseline/`             | ctrl  | `damage_gradient`, `n_d = grad(d)/\|grad(d)\|`       | damage, bounded `[0.008, 0.999]`         |
| `regularized_normal/`   | A     | `damage_gradient` **+ regularization** `ε = 1e-8`   | damage (= baseline)                      |
| `strained-based/`       | A     | `principal_strain`, `n_F = e₁` (max principal)      | damage (= baseline)                      |
| `porosity_bounded/`     | B     | `damage_gradient` (= baseline)                      | damage, **upper bound `0.065`**          |
| `porosity-strain-based/`| B     | `damage_gradient` (= baseline)                      | **`strain`**, `φ = φ₀ + ε₁`              |

References: damage-gradient normal — Heider 2021 eq. 46; strain-based normal —
Liu et al. 2024 *CMAME* **429**:117165 eqs. 29–30; strain-based porosity —
Liu et al. 2024 eq. 40.

## `[elasticity]` block differences (everything else identical)
- `baseline/`: `crack_normal_source = damage_gradient` (no regularization);
  `[porosity_damaged]` default `damage` law, `porosity_upper_bound = 0.999`.
- `regularized_normal/`: adds `regularize_crack_normal = true` +
  `crack_normal_regularization = 1e-8`.
- `strained-based/`: `crack_normal_source = principal_strain`.
- `porosity_bounded/`: `porosity_upper_bound = 0.065` (caps the damage law).
- `porosity-strain-based/`: `porosity_update_model = strain` +
  `strain_property = mechanical_strain`, `porosity_upper_bound = 0.999`.

Verify the minimal diffs:
```bash
diff baseline/elasticity_E1d25.i strained-based/elasticity_E1d25.i      # only crack_normal_source
diff baseline/elasticity_E1d25.i regularized_normal/elasticity_E1d25.i  # only +regularize keys
diff baseline/elasticity_E1d25.i porosity_bounded/elasticity_E1d25.i    # only porosity_upper_bound
diff baseline/elasticity_E1d25.i porosity-strain-based/elasticity_E1d25.i # only porosity keys (+ header note)
```

## Two meshes per case (resolution / convergence)
Each case ships two problems that differ **only** in mesh resolution and the
matching static initial condition:

| Problem (per case)                | Mesh                          | Static IC (SolutionUserObject)      |
|-----------------------------------|-------------------------------|-------------------------------------|
| `elasticity_E1d25.i`              | `../../../2d_mesh/2d_mesh.msh` (lc = 2e-4) | `../static_solve_out.e`        |
| `elasticity_E1d25_refined.i`      | `../../../2d_mesh/2d_mesh_refined.msh` (lc = 5e-5, ≈16× elements) | `../static_solve_refined_out.e` |

Each main app drives its own fracture sub-app (`fracture_E1d25.i` /
`fracture_E1d25_refined.i`) on the same mesh.

## Shared static initial conditions
All cases read the **same two** static-solve solutions from the parent
`permeability_formula/` level (no per-case copies):
- `static_solve_out.e`         — base mesh, from `static_solve.i`
- `static_solve_refined_out.e` — refined mesh, from `static_solve_refined.i`

Both are generated once (see *Running* below) and reused by all five cases via
`SolutionUserObject` (`disp_x disp_y pp elastic_strain_*`). Because each refined
problem reads a static IC computed **on the same refined mesh**, the nodal
transfer is exact (no interpolation).

## Consistent diagnostics (all 10 runs)
Every elasticity input writes the same Exodus `show` fields so the cases are
directly comparable:
- `bulk_modulus_degraded_aux` — `K = (1/9) I:C:I` from the degraded SPECTRAL
  elastic tangent (≈ `K` where `d→0`/compressed, drops in damaged tension).
- `stress_00 stress_01 stress_11` — Cauchy effective-stress components.
- `effective_perm00_aux effective_perm11_aux effective_perm01_aux` — fracture
  permeability tensor (`MONOMIAL CONSTANT`).
- plus `d`, velocities, `pp`, `psie_active_enhanced`, `biot_modulus_aux`,
  `biot_coefficient_aux`, `porosity_aux`.

Run controls (identical across all cases): `permeability_model = normal_strain`,
`decomposition = SPECTRAL` (required for `effective_perm` to update),
`dt = 1e-8`, `NewmarkBeta`, **`end_time = 10e-5` (10 pulses)**.

## Expected distinguishing effects
**Group A** — at the fully-damaged crack core (`d → 1`, `grad d → 0`):
- `strained-based/`: `n_F = e₁` stays a unit vector ⇒ projector `(I − n_F⊗n_F)`
  keeps fracture permeability **anisotropic** (normal direction suppressed).
- `regularized_normal/`: `n_d → 0` ⇒ `(I − n_d⊗n_d) → I` ⇒ permeability becomes
  **isotropic** at the core.
- `baseline/`: unregularized `damage_gradient` — `n_d` ill-conditioned where
  `grad d → 0` (the reason regularization / strain normals are studied).

Compare `effective_perm*` (and the resulting `pp`/aperture/leak-off).

**Group B** — porosity response:
- `porosity_bounded/`: damage-driven `φ` saturates at the `0.065` cap.
- `porosity-strain-based/`: `φ = φ₀ + ε₁` tracks the opening strain (reversible),
  capped only at `0.999`. Compare `porosity_aux` and porosity-dependent outputs
  (Biot modulus, fluid energy). Note: the `*_static` energy constants are
  calibrated for the baseline porosity model and should be re-derived before
  using them to interpret the strain-porosity energy balance (the
  `porosity_aux` field comparison itself is unaffected).

## Running

### 1. Static initial conditions (once)
- **Base mesh:** `static_solve_out.e` already exists in this directory.
- **Refined mesh:** generate `static_solve_refined_out.e` from
  `static_solve_refined.i` (the mesh-refined twin of `static_solve.i`). It can be
  produced **locally** (a single Steady solve) and the resulting Exodus file is
  portable, so no separate cluster static job is required:
  ```bash
  # from permeability_formula/, in the moose conda env
  ../../../../farms-opt -i static_solve_refined.i        # writes static_solve_refined_out.e
  # (or: mpiexec -n <N> ../../../../farms-opt -i static_solve_refined.i)
  ```
  `submit_static_refined.sbatch` remains as a cluster fallback if local
  generation is not possible.

### 2. Elasticity runs
Locally, run any case directly, e.g.:
```bash
farms-opt -i baseline/elasticity_E1d25.i            # base mesh
farms-opt -i baseline/elasticity_E1d25_refined.i    # refined mesh
```
On the cluster, submit everything with the batch driver one level up:
```bash
cd ..                       # 2d_hydromech/
./submit_all.sh --dry-run   # list the discovered jobs (includes both meshes)
./submit_all.sh             # submit (static ICs must already be in place)
```
`submit_all.sh` discovers both `submit_elasticity.sbatch` and
`submit_elasticity_refined.sbatch` for every case and submits each from its case
directory (logs/outputs land next to the inputs). The static solves are **not**
submitted by the driver — upload the two `static_solve*_out.e` files into this
directory first.

> The sbatch scripts use absolute `/scratch2/...` paths and per-case job names;
> update the path prefix to your checkout location before submitting.
