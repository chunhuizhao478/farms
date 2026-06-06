# Expected values — `damaged_porosity` material tests

Hand-computed gold values for the `ElkPorousFlowDamagedPorosity` unit tests. Every test
imposes an **element-uniform** field (piecewise-constant damage, or piecewise-linear
displacement with element-constant slope) so the porosity is exact per element with no
quadrature averaging. Constants: `phi_0 = 0.008`.

The clamp applied to both laws is `phi = clamp(phi_raw, lower, upper)`.

---

## `damage_model.i` — default `damage` law

Law: `phi_raw = phi_0 + (1 - phi_0) * (1 - (1 - d)^2)`, with `lower = 0.008`, `upper = 0.999`.

Damage is piecewise-constant with breakpoints on the element boundaries (`x = 1/3, 2/3`),
so each element carries one exact `d`:

| element | x-range     | d     | (1-d)^2                | 1-(1-d)^2   | phi = 0.008 + 0.992·(1-(1-d)^2) |
|--------:|-------------|-------|------------------------|-------------|--------------------------------|
| 0       | [0, 1/3]    | 1/6   | (5/6)^2 = 0.694444444  | 0.305555556 | **0.311116622**                |
| 1       | [1/3, 2/3]  | 1/2   | 0.25                   | 0.75        | **0.752000000**                |
| 2       | [2/3, 1]    | 5/6   | (1/6)^2 = 0.027777778  | 0.972222222 | **0.972443333**                |

Sample points (element centers): `phi_elem0 @ (0.1666667, 0.05)`,
`phi_elem1 @ (0.5, 0.05)`, `phi_elem2 @ (0.8333333, 0.05)`.
All three are within `[0.008, 0.999]`, so no clamping occurs.

---

## `strain_model.i` — `strain` law (Liu et al. 2024 eq. 40)

Law: `phi_raw = phi_0 + eps_1`, with `eps_1` the maximum principal strain of
`mechanical_strain`, `lower = 0.008`, `upper = 0.999`.

`disp_x` is continuous piecewise-linear with nodal values `[0, 1e-3/3, ...]` chosen so the
element gradient `eps_xx` equals the target slope; `disp_y = 0`. In 2D plane strain the
strain tensor per element is `diag(eps_xx, 0, 0)`, whose eigenvalues (ascending) are
`[min(eps_xx,0), 0, max(eps_xx,0)]`, so `eps_1 = max(eps_xx, 0)`.

| element | slope eps_xx | eigenvalues        | eps_1 = max(eps_xx,0) | phi_raw = 0.008 + eps_1 | phi (clamped)   |
|--------:|--------------|--------------------|-----------------------|-------------------------|-----------------|
| 0       | 1e-3         | [0, 0, 1e-3]       | 1e-3                  | 0.009                   | **0.009**       |
| 1       | 2.0          | [0, 0, 2.0]        | 2.0                   | 2.008                   | **0.999** (upper) |
| 2       | -0.5         | [-0.5, 0, 0]       | 0                    | 0.008                   | **0.008**       |

Sample points: `phi_tension @ (0.1666667, 0.05)`,
`phi_clamp_upper @ (0.5, 0.05)`, `phi_compression @ (0.8333333, 0.05)`.

Notes:
- Element 1 demonstrates the **upper clamp** (2.008 → 0.999).
- Element 2 demonstrates that **pure compression yields `phi = phi_0`** — the most-tensile
  principal strain is 0 (`eps_zz = eps_yy = 0`), not the negative `eps_xx`. This is the
  physically intended behavior of eq. (40), not the lower clamp.
- `rel_err = 1e-6` in the test absorbs the ~1e-12 perturbation from finite-digit breakpoints.

### `strain_lower_clamp` (reuses `strain_model.i`, `porosity_lower_bound = 0.05`)

With `lower = 0.05`:

| element | phi_raw | phi (clamped to [0.05, 0.999]) |
|--------:|---------|--------------------------------|
| 0       | 0.009   | **0.05** (lower clamp)         |
| 1       | 2.008   | **0.999** (upper clamp)        |
| 2       | 0.008   | **0.05** (lower clamp)         |

This exercises the **lower clamp**: in 2D the raw strain porosity can never fall below
`phi_0`, so a lower bound above `phi_0` is the way to drive the lower-clamp branch.

---

## `strain_missing_strain_error` — error path

`strain_missing_strain.i` selects `porosity_update_model = strain` but supplies no
`ComputeSmallStrain`, so `mechanical_strain` is undefined. The run must abort with an error
mentioning `mechanical_strain` (MOOSE "material property … not defined"), confirming the
strain law fails loudly instead of defaulting to a zero strain.

---

## Regenerating gold

Gold CSVs are produced by running the inputs with a working `farms-opt` (see `regold.sh`).
After regolding, confirm the numbers above appear in `gold/*.csv`. The local macOS build
links `libfarms-opt.dylib` but the final `farms-opt` executable does not link in the conda
moose env on macOS 26.x; generate gold and run `./run_tests` on a platform where the
executable links (e.g. CI / Linux).
