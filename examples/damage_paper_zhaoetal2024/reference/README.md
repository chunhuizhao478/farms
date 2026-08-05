# Reference fingerprint

`fingerprint_alpha_B.csv` is the pass/fail target for this example. A run reproduces the
published case when it matches this table — not merely when it completes.

## Where these numbers come from

They were measured directly from the Exodus output of the original published run:

| | |
| --- | --- |
| file | `test_planarfault_main_out.e` |
| size | 11,940,860,476 bytes |
| run date | 2023-11-20 19:44:55 |
| MOOSE | git commit `9de1ffd476` (2023-11-09) |
| PETSc | 3.20.0 |
| farms | branch `cdbm`, commit `99b41e9c` |
| frames | 31, at t = 0.0, 0.1, ... 3.0 s |

The provenance above was read out of the file's own `info_records` NetCDF variable, which MOOSE
writes into every Exodus output. It records the full command line, the framework versions, the
run timestamp, and the complete resolved input deck including defaults the deck never set. If you
ever need to identify which code produced an orphaned `.e` file, read that variable first — do
not infer it from directory names.

## What the columns mean

`alpha` and `B` are the CDBM damage and breakage variables. Both are `CONSTANT MONOMIAL`, i.e.
Exodus **element** variables, and the statistics pool both element blocks (372,298 + 372,584
triangles). `n_gt_1em6` counts elements with value greater than 1e-6.

## Tolerance policy

- **`n_gt_1em6` must match exactly.** These are integers produced by threshold crossings, so any
  drift in the fields moves them. They are the sharpest discriminator available and the reason
  this fingerprint is meaningful rather than decorative.
- **`max` and `mean` must agree to 6 significant figures.** They are stored here at 6 significant
  figures, which is the precision at which they were measured.
- `min` is 0 everywhere and is retained only as a sanity check that the field is non-negative.

## The two signatures that matter most

1. **`alpha` reaches exactly 1.0 from t = 1.7 s onward.** A build carrying the wrong critical-damage
   calibration caps `alpha` at **0.952169427** instead. That single number is the fastest way to
   detect a miscalibrated build — see `../README.md`, "Known pitfalls".
2. **`B` is identically zero through t = 1.3 s, then first appears at t = 1.4 s in exactly 15
   elements with max 0.605578.** Fifteen elements at a specific value is a far more discriminating
   check than any bulk field norm.

## What this fingerprint does NOT prove

Matching output **file size** proves nothing about the physics. A rerun of this case with a
differently-calibrated build produced an Exodus of 11,940,877,188 bytes against the reference's
11,940,860,476 — a difference of 16,712 bytes in 11.94 GB, or 0.00014% — while containing
materially different damage (45,179 broken elements at t = 3.0 s instead of 48,643). File size
tracks mesh, variable set, and frame count. It says nothing about the values inside.

Use `../verify_fingerprint.py`. Do not eyeball file sizes.
