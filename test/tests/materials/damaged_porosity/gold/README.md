# gold/ — placeholder

The CSV gold files referenced by `../tests` are produced by running the inputs with a
working `farms-opt` executable. They are **not** checked in yet because the local macOS 26.x
build links `libfarms-opt.dylib` but not the final executable (see
`../../../../memory`/toolchain note and `../EXPECTED_VALUES.md`).

To populate this directory on a platform where the executable links:

```
../regold.sh           # or: FARMS_APP=/path/to/farms-opt ../regold.sh
```

Expected files after regolding:
- `damage_model_out.csv`         — phi = 0.311116622, 0.752, 0.972443333
- `strain_model_out.csv`         — phi = 0.009, 0.999, 0.008
- `strain_lower_clamp_out.csv`   — phi = 0.05, 0.999, 0.05

Verify the values against `../EXPECTED_VALUES.md` before committing the CSVs (and delete
this README once the gold is in place).
