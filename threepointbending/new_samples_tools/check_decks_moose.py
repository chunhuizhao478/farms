#!/usr/bin/env python3
"""
check_decks_moose.py -- run the new-sample decks through the local farms-opt
without solving (MOOSE 94b6528, the commit of Frontera's moose94b6528 build).

WHAT
    1. `--check-input` on every code_<family>/case_*/elasticity.i: MOOSE reads the
       .msh, runs the deck's mesh generators (elastic pads under the supports and
       the load, support/loading nodesets), builds every object and stops before
       the first time step.  Proves each deck and its mesh work together.
    2. `--split-mesh 400` (the serial first step of every Frontera job) on the
       largest mesh, recording wall time and peak memory.

STEPS
    run_moose()   farms-opt under mesh_new_samples.run_guarded() (RSS / free-memory
                  kill-switch), working directory = a scratch folder, so JIT caches
                  and split files never land in the repository

CHECKS
    exit status 0 and no "*** ERROR" line in the log.

OUTPUT
    logs/check_decks_moose.tsv, logs/moose_<case>_<mode>.log

USAGE
    python check_decks_moose.py --scratch DIR [--only SUB1,SUB2] [--split-only | --check-only]
"""

import argparse
import os
import shutil
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import make_new_sample_cases as cases  # noqa: E402
import make_new_sample_geos as geos  # noqa: E402
import mesh_new_samples as mesher  # noqa: E402

FARMS_OPT = "/Users/chunhuizhao/projects/farms_cdms/farms-opt"   # Mac build, MOOSE 94b6528
MOOSE_ENV_BIN = "/Users/chunhuizhao/miniforge/envs/moose/bin"
SDKROOT = "/Library/Developer/CommandLineTools/SDKs/MacOSX15.4.sdk"  # JIT links fail with the 26.5 SDK
N_SPLITS = 400                                                    # as NUM_SPLITS in the jobs
SUMMARY_TSV = os.path.join(mesher.LOG_DIR, "check_decks_moose.tsv")


def moose_env():
    env = dict(os.environ)
    env["PATH"] = MOOSE_ENV_BIN + os.pathsep + env.get("PATH", "")
    env["SDKROOT"] = SDKROOT
    return env


def run_moose(case, mode, scratch):
    """Run one check ('check' or 'split') for one case; return (ok, status, wall, rss, free)."""
    deck = os.path.join(geos.TPB_DIR, f"code_{case['family']}", cases.case_dir_name(case), "elasticity.i")
    work = os.path.join(scratch, f"{geos.case_tag(case)}_{mode}")
    os.makedirs(work, exist_ok=True)
    cmd = [FARMS_OPT, "-i", deck, "--color", "off"]
    if mode == "check":
        cmd += ["--check-input"]
    else:
        cmd += ["--split-mesh", str(N_SPLITS), "--split-file", os.path.join(work, "foo.cpr")]
    log_path = os.path.join(mesher.LOG_DIR, f"moose_{geos.case_tag(case)}_{mode}.log")
    status, wall, rss, free = mesher.run_guarded(cmd, log_path, cwd=work, env=moose_env())
    with open(log_path) as fh:
        errors = [line for line in fh if "*** ERROR" in line]
    ok = status == "ok" and not errors
    if errors and status == "ok":
        status = "MOOSE error"
    shutil.rmtree(work, ignore_errors=True)
    return ok, status, wall, rss, free


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scratch", required=True, help="scratch folder for working directories and split files")
    ap.add_argument("--only", default="", help="comma-separated substrings; a case matching any is run")
    ap.add_argument("--split-only", action="store_true", help="only the split of the largest selected mesh")
    ap.add_argument("--check-only", action="store_true", help="skip the split")
    args = ap.parse_args()

    patterns = args.only.split(",")
    selected = [c for c in geos.CASES if any(p in geos.mesh_basename(c) for p in patterns)]
    jobs = [] if args.split_only else [(c, "check") for c in selected]
    if not args.check_only:
        largest = max(selected, key=lambda c: os.path.getsize(geos.geo_path(c)[:-4] + ".msh"))
        jobs.append((largest, "split"))

    os.makedirs(mesher.LOG_DIR, exist_ok=True)
    failures = 0
    with open(SUMMARY_TSV, "a") as summary:
        for case, mode in jobs:
            ok, status, wall, rss, free = run_moose(case, mode, args.scratch)
            failures += not ok
            row = f"{geos.case_tag(case)}\t{mode}\t{status}\t{wall:.0f}\t{rss:.2f}\t{free:.0f}"
            summary.write(row + "\n")
            summary.flush()
            print(f"{'OK ' if ok else 'BAD'} {mode:5s} {geos.case_tag(case):45s} {status}: {wall:.0f} s, "
                  f"peak RSS {rss:.2f} GB, min free {free:.0f} %", flush=True)
    print("DONE" if failures == 0 else f"DONE with {failures} failures", flush=True)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
