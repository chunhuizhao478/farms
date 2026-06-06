# Crack-normal comparison: strain-based vs. regularized damage-gradient

These cases isolate the **crack-normal definition** `n_F` in the `normal_strain`
permeability model. Everything else in a given case (mesh, BCs, loading,
`perm_exponent`, `characteristic_length_type = element_size`, aperture
`w_c = h_c·|1+ε_nn|`, threshold, `d^b` weighting) is unchanged — only
`crack_normal_source` differs.

| Case dir            | Base case copied from        | Crack normal `n_F`                              | Reference                          |
|---------------------|------------------------------|-------------------------------------------------|------------------------------------|
| `strained-based/`   | `undrained/em_0p005/`        | `n_F = e₁`, eigenvector of max principal strain | Liu et al. 2024 CMAME eqs. 29–30   |
| `regularized_normal/`| top-level `permeability_formula/` | `n_d = grad(d)/(|grad(d)| + ε)`, `ε = 1e-8` | Heider 2021 eq. 46 + regularization |

`[elasticity]` block difference:
- `strained-based/elasticity_E1d25.i`: `crack_normal_source = principal_strain`.
- `regularized_normal/elasticity_E1d25.i`: `crack_normal_source = damage_gradient`
  + `regularize_crack_normal = true` + `crack_normal_regularization = 1e-8`.

## ⚠ Base-case caveat (for a strict A/B)
`strained-based/` is built from the `undrained/em_0p005/` case, while
`regularized_normal/` is built from the top-level `permeability_formula/` case —
**different base setups**. They are not a same-base A/B pair. For an apples-to-apples
comparison of the normal definition at em_0p005 resolution, compare `strained-based/`
against an em_0p005 baseline:
- the existing `undrained/em_0p005/` (uses `crack_normal_source = damage_gradient`,
  no regularization), or
- a regularized em_0p005 variant (`crack_normal_source = damage_gradient` +
  `regularize_crack_normal = true`) — handled separately in another worktree.

## Expected distinguishing effect
At the fully-damaged crack core (`d → 1`, `grad(d) → 0`):
- **strain-based** (`strained-based/`): `n_F = e₁` stays a well-defined unit vector,
  so the tangential projector `(I − n_F⊗n_F)` keeps the fracture permeability
  **anisotropic** (suppressed in the crack-normal direction).
- **regularized damage-gradient** (`regularized_normal/`): `n_d → 0`, so
  `(I − n_d⊗n_d) → I` and the fracture permeability becomes **isotropic** at the core.

Compare the `effective_perm` field (and the resulting pressure/aperture/leak-off)
between runs.

## Running (cluster)
Each case reuses the shared `../static_solve_out.e` (initial conditions via the
`SolutionUserObject` in `elasticity_E1d25.i`) and the shared mesh
`../../../2d_mesh/2d_mesh.msh`. Generate `static_solve_out.e` once at the
`permeability_formula/` level (run `static_solve.i` / `submit_static.sbatch`), then:

```bash
sbatch strained-based/submit_elasticity.sbatch
sbatch regularized_normal/submit_elasticity.sbatch
```

> Note: `permeability_model = normal_strain` only updates `effective_perm` under
> `decomposition = SPECTRAL` (already set in these inputs). The submit scripts use
> absolute `/scratch2/...` paths and a per-case job name; update the path prefix to
> your checkout location before submitting.
