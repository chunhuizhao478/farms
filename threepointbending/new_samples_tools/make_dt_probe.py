#!/usr/bin/env python3
"""
make_dt_probe.py -- a coarse-dt probe of one new-sample case, to find when the
solve becomes hard (crack nucleation), i.e. where a small dt is really needed.

WHAT
    Copies code_<family>/case_<family>_3d_<tag>/ into a sibling folder
    case_<family>_3d_<tag>_dt0p1_end50/ (production outputs and checkpoints are
    never touched) and changes only:
      - dt_fn       0.1 for all t (production: 0.1 to t = 14, then 0.001)
      - end_time    50 (production: 20)
      - sync_times  every 0.5 up to 50 (Exodus frames)
    and adds per-step CSV postprocessors: max d, fixed-point / nonlinear / linear
    iterations of the step, and elapsed wall time.  Solver settings are unchanged,
    so "hard" means the same thing as in production.
    Also writes jobs_<family>_dtprobe/job_<tag>_dt0p1_end50.sbatch: development
    queue, PROBE_NODES nodes x 40 ranks (production: 10 x 40), 2 h, own mesh split.

STEPS
    1. probe_deck_texts()  production deck + asserted substitutions
    2. probe_job_text()    production job + development queue, probe case path
    3. check_with_moose()  --check-input with the local farms-opt (run_guarded)

CHECKS
    Every substitution must match exactly once; --check-input must exit 0.

USAGE
    python make_dt_probe.py TAG [--check]
      TAG e.g. xL3.00mm_xR2.00mm_y3.6mm_radi1.0mm (a case_tag of CASES)
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import make_new_sample_cases as cases  # noqa: E402
import make_new_sample_geos as geos  # noqa: E402

PROBE_DT = 0.1
PROBE_END_TIME = 50
PROBE_SYNC_STEP = 0.5
PROBE_SUFFIX = "_dt0p1_end50"
PROBE_NODES = 20           # development queue (max 40); production uses 10
PROBE_RANKS = 40 * PROBE_NODES
PROBE_POSTPROCESSORS = """
# --- dt probe diagnostics (one CSV row per time step) ---
[Postprocessors]
  [max_d]
    type = NodalExtremeValue
    variable = d
  []
  [fixed_point_its]
    type = NumFixedPointIterations
  []
  [nonlinear_its]
    type = NumNonlinearIterations
  []
  [linear_its]
    type = NumLinearIterations
  []
  [wall_time]
    type = PerfGraphData
    section_name = Root
    data_type = TOTAL
  []
[]
"""


def find_case(tag):
    for case in geos.CASES:
        if geos.case_tag(case) == tag:
            return case
    raise SystemExit(f"unknown case tag {tag!r}")


def probe_deck_texts(case):
    """(elasticity.i, fracture.i) of the probe."""
    elastic, fracture = cases.deck_texts(case)
    elastic = cases.substitute(elastic, "x = '0   14'\n    y = '0.1 0.001'",
                               f"x = '0   {PROBE_END_TIME}'\n    y = '{PROBE_DT} {PROBE_DT}'", "probe dt_fn")
    elastic = cases.substitute(elastic, "  end_time = 20\n", f"  end_time = {PROBE_END_TIME}\n", "probe end_time")
    start = elastic.index("sync_times = '")
    end = elastic.index("'", start + len("sync_times = '"))
    n = int(round(PROBE_END_TIME / PROBE_SYNC_STEP))
    times = " ".join(f"{PROBE_SYNC_STEP * k:g}" for k in range(1, n + 1))
    elastic = elastic[:start] + f"sync_times = '{times}" + elastic[end:]
    elastic = cases.substitute(elastic, "[Outputs]\n", PROBE_POSTPROCESSORS + "\n[Outputs]\n  csv = true\n", "probe outputs")
    return elastic, fracture


def probe_job_text(case):
    tag, fam = geos.case_tag(case), case["family"]
    t = cases.job_text(case)
    t = t.replace(f"3dptrbend_{fam}_{tag}", f"3dptrbend_{fam}_{tag}{PROBE_SUFFIX}")
    t = cases.substitute(t, "#SBATCH -p normal ", "#SBATCH -p development ", "probe queue")
    t = cases.substitute(t, "#SBATCH -t 48:00:00 ", "#SBATCH -t 02:00:00 ", "probe wall time")
    t = cases.substitute(t, "#SBATCH -N 10 ", f"#SBATCH -N {PROBE_NODES} ", "probe nodes")
    t = cases.substitute(t, "#SBATCH -n 400 ", f"#SBATCH -n {PROBE_RANKS} ", "probe ranks")
    t = cases.substitute(t, "NUM_SPLITS=400", f"NUM_SPLITS={PROBE_RANKS}", "probe split count")
    t = cases.substitute(t, f'CASE_NAME="{tag}"', f'CASE_NAME="{tag}{PROBE_SUFFIX}"', "probe case name")
    return t


def check_with_moose(deck, scratch):
    """--check-input with the local farms-opt under the memory guard; returns the status."""
    import check_decks_moose as checker
    import mesh_new_samples as mesher
    os.makedirs(scratch, exist_ok=True)
    log = os.path.join(mesher.LOG_DIR, "moose_dtprobe_check.log")
    status, wall, rss, free = mesher.run_guarded([checker.FARMS_OPT, "-i", deck, "--color", "off", "--check-input"],
                                                 log, cwd=scratch, env=checker.moose_env())
    with open(log) as fh:
        if any("*** ERROR" in line for line in fh):
            status = "MOOSE error"
    print(f"--check-input: {status}, {wall:.0f} s, peak RSS {rss:.2f} GB (log {log})")
    return status


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("tag")
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--scratch", default="/tmp/dtprobe_check")
    args = ap.parse_args()
    case = find_case(args.tag)
    fam = case["family"]
    case_dir = os.path.join(geos.TPB_DIR, f"code_{fam}", cases.case_dir_name(case) + PROBE_SUFFIX)
    elastic, fracture = probe_deck_texts(case)
    cases.write(os.path.join(case_dir, "elasticity.i"), elastic)
    cases.write(os.path.join(case_dir, "fracture.i"), fracture)
    job = os.path.join(geos.TPB_DIR, f"jobs_{fam}_dtprobe", f"job_{args.tag}{PROBE_SUFFIX}.sbatch")
    cases.write(job, probe_job_text(case))
    print(f"wrote {os.path.relpath(case_dir, geos.TPB_DIR)}/ and {os.path.relpath(job, geos.TPB_DIR)}")
    if args.check and check_with_moose(os.path.join(case_dir, "elasticity.i"), args.scratch) != "ok":
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
