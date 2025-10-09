# Repository Guidelines

## Project Structure & Module Organization
- `src/` holds the primary MOOSE application code (AuxKernels, Kernels, Materials, BCs). Mirror types in `include/` for headers; keep new components in matching subdirectories.
- `test/` mirrors the production layout (`test/src`, `test/include`, `test/tests`) and supplies harness fixtures and `.i` decks; place new regression decks in an existing physics folder or create a descriptive directory.
- `examples/` contains published scenarios (e.g., `examples/benchmark_tpv2052D`); use these as templates for shareable input files, and document any new case-specific assets.
- `scripts/` hosts helper utilities; keep prototypes under version control only when they are reusable.
- The bundled `moose/` submodule is the default framework; set `MOOSE_DIR` if pointing to an external checkout.

## Build, Test, and Development Commands
- `make -j8` builds the optimized binary `farms-opt`; use `METHOD=dbg make -j4` when you need debug symbols (`farms-dbg`).
- `./farms-opt -i examples/benchmark_tpv2052D/tpv2052D.i` runs a representative case; prefer relative paths rooted at the repository.
- `./run_tests --re aftershock2D` executes the MOOSE TestHarness with a regex filter; omit `--re` to run the full regression suite and combine with `-j 4` to parallelize.

## Coding Style & Naming Conventions
- C++ code follows the repository `.clang-format` (LLVM base, 2-space indentation, Allman braces, 100-character columns); run `clang-format` before committing non-trivial changes.
- Class names use PascalCase (`ComputeDamageBreakageStress3DSlipWeakening`), method overrides stay camelCase, and member data uses the `m_` prefix consistent with MOOSE.
- Input files use lowercase with underscores (`aftershock2d.i`); keep directory names descriptive and short.

## Testing Guidelines
- Tests extend the MOOSE harness; add `.i` decks and gold files in `test/tests/<feature>/`, naming them after the behavior under evaluation (`damage_breakage_stability.i`).
- Cover new kernels or materials with steady and dynamic runs when applicable, and document non-standard parameters in comments for future reproducibility.
- Use `./run_tests --failed-tests-report` before opening a pull request to confirm the suite and capture a concise summary for reviewers.

## Commit & Pull Request Guidelines
- Follow the observed `feature_scope: imperative summary` pattern (`aftershock2D: adjust dashpot BC strength`); keep subject lines under 72 characters and explain deltas in the body when needed.
- Group related changes per commit (code + deck + docs together).
- Pull requests should outline the physics motivation, list affected input files, reference GitHub issues if available, and attach relevant plots or log excerpts from `test/tests`.
- Verify that the branch rebases cleanly on `main` and that CI (or local `run_tests`) passes before requesting review.

## Environment & Configuration Tips
- Ensure PETSc and MOOSE prerequisites are installed; source the MOOSE environment script before building if you rely on a centralized `MOOSE_DIR`.
- Capture complex solver settings in `.yaml` or `.i` files committed under `examples/` or `doc/`, and keep paths portable by using environment variables where needed.
