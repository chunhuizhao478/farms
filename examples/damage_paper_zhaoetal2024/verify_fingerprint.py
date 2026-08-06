#!/usr/bin/env python3
"""Check a run of this example against the published fingerprint.

    python3 verify_fingerprint.py test_planarfault_main_out.e

Exit 0 only if every frame matches. "It ran" is not the success criterion; matching
these numbers is.

Reads the Exodus frame by frame. That is deliberate: the output is ~12 GB, an mmap of
that size is refused under a typical login-node address-space limit, and
scipy.io.netcdf_file without mmap loads the entire record array into RAM. Run this as a
small batch job if your login node is tight.

Do NOT run this against an output file while the job producing it is still writing:
on a parallel filesystem the reader forces write-lock revocation and can stall the
writer for many minutes.
"""
import argparse
import csv
import pathlib
import sys

import numpy as np
from scipy.io import netcdf_file

# alpha and B are CONSTANT MONOMIAL, i.e. Exodus *element* variables.
VARMAP = {"alpha": "alpha_in", "B": "B_in"}
SIGFIGS = 6


def sig(x, n=SIGFIGS):
    """Round to n significant figures, so comparisons are at the precision recorded."""
    if x == 0 or not np.isfinite(x):
        return 0.0
    return float(f"%.{n}g" % x)


def decode_names(var):
    out = []
    for row in var[:]:
        s = b"".join(bytes([c]) if isinstance(c, (int, np.integer)) else bytes(c)
                     for c in row if c not in (0, b"\x00"))
        out.append(s.decode(errors="replace").strip())
    return out


def stats_for_frame(f, var_index, nblk, k):
    vmax, vmin, tot, cnt, nz = -np.inf, np.inf, 0.0, 0, 0
    for b in range(1, nblk + 1):
        key = f"vals_elem_var{var_index}eb{b}"
        if key not in f.variables:
            continue
        a = np.asarray(f.variables[key][k, :], dtype=float)
        if a.size == 0:
            continue
        vmax = max(vmax, float(a.max()))
        vmin = min(vmin, float(a.min()))
        tot += float(a.sum())
        cnt += a.size
        nz += int((a > 1e-6).sum())
    return vmax, vmin, (tot / cnt if cnt else 0.0), nz


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("exodus", help="output .e from this example")
    ap.add_argument("--fingerprint",
                    default=str(pathlib.Path(__file__).parent / "reference" / "fingerprint_alpha_B.csv"))
    args = ap.parse_args()

    expected = {}
    with open(args.fingerprint) as fh:
        for row in csv.DictReader(fh):
            expected[(round(float(row["time"]), 4), row["var"])] = row

    f = netcdf_file(args.exodus, "r", mmap=True)
    names = decode_names(f.variables["name_elem_var"])
    times = np.asarray(f.variables["time_whole"][:], dtype=float)
    nblk = f.dimensions["num_el_blk"]

    print(f"file   : {args.exodus}")
    print(f"frames : {len(times)}  t=[{times[0]:g} .. {times[-1]:g}]")

    failures = []
    if len(times) != 31:
        failures.append(f"expected 31 frames, found {len(times)} "
                        f"(a short run means the deck's num_steps was left active)")

    print(f"\n{'time':>6} {'var':>6} {'max':>13} {'expected':>13} {'n':>8} {'expected':>8}  status")
    for label, exo_name in VARMAP.items():
        if exo_name not in names:
            failures.append(f"element variable {exo_name} not present in the output")
            continue
        vi = names.index(exo_name) + 1
        for k, t in enumerate(times):
            key = (round(float(t), 4), label)
            if key not in expected:
                key = (round(round(float(t), 1), 4), label)   # tolerate 0.30000000000004
            if key not in expected:
                failures.append(f"t={t:g} {label}: no fingerprint row")
                continue
            exp = expected[key]
            vmax, _vmin, mean, nz = stats_for_frame(f, vi, nblk, k)
            e_max, e_mean, e_n = float(exp["max"]), float(exp["mean"]), int(exp["n_gt_1em6"])

            bad = []
            if nz != e_n:
                bad.append(f"count {nz} != {e_n}")
            if sig(vmax) != sig(e_max):
                bad.append(f"max {sig(vmax)} != {sig(e_max)}")
            if sig(mean) != sig(e_mean):
                bad.append(f"mean {sig(mean)} != {sig(e_mean)}")

            status = "ok" if not bad else "FAIL: " + "; ".join(bad)
            if bad:
                failures.append(f"t={t:.1f} {label}: " + "; ".join(bad))
            print(f"{t:6.1f} {label:>6} {vmax:13.6g} {e_max:13.6g} {nz:8d} {e_n:8d}  {status}")

    f.close()

    print()
    if not failures:
        print("PASS - reproduces the published fingerprint.")
        return 0

    print(f"FAIL - {len(failures)} mismatch(es):")
    for x in failures[:20]:
        print("  -", x)
    if len(failures) > 20:
        print(f"  ... and {len(failures) - 20} more")

    # The single most likely cause, called out by name.
    alpha_capped = any("max 0.952169" in x for x in failures)
    if alpha_capped:
        print("\nDIAGNOSIS: alpha is capping at 0.952169427 instead of reaching 1.0.")
        print("That is the signature of the wrong alpha_cr calibration -- the breakage")
        print("kernels are using the 2.73e9 default while this case needs 3.204e10.")
        print("Check that both BreakageVarForcingFuncDevOld blocks in")
        print("test_planarfault_sub.i set: alphacr_calibration = stiff_3p204e10")
    return 1


if __name__ == "__main__":
    sys.exit(main())
