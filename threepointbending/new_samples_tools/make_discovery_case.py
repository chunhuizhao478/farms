#!/usr/bin/env python3
"""
make_discovery_case.py -- a fresh-start deck for one new-sample case on USC Discovery.

WHAT
    Discovery starts each case from t = 0, so it can skip the 6,000 fine steps of the Frontera
    production runs between t = 14 and 20: the dt probes showed dt = 0.1 is accurate until the
    crack (probe max d within 2-3 % of production).  Writes
    code_<family>/case_<family>_3d_<tag>_discovery/{elasticity.i, fracture.i} with
      - dt_fn       0.1 from t = 0 to t_fine, then 0.001  (t_fine as make_restart_decks.py:
                    crack - 1.5, down to 0.5, >= 20)
      - end_time    50 (as the Frontera restarts; the wall time ends the run first)
      - sync_times  every 1.0 up to t_fine (and at t_fine), every 0.1 after it
      - per-step CSV: max d, fixed-point / nonlinear / linear iterations, wall time (PerfGraphData)
        and peak process memory (MemoryUsage), to measure the step rate and memory need
    Physics and solver settings are the production deck's (code_2holes template).  Hourly
    checkpoints stay on (MOOSE default), so a later job can continue with --recover.

STEPS
    1. discovery_deck_texts()  production deck + asserted substitutions
    2. write the case folder

CHECKS
    Every substitution must match exactly once; check with check_decks_moose-style --check-input.

USAGE
    python make_discovery_case.py TAG
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import make_new_sample_cases as cases  # noqa: E402
import make_new_sample_geos as geos  # noqa: E402
import make_restart_decks as restart  # noqa: E402

SUFFIX = "_discovery"
DIAGNOSTICS = """
# --- per-step diagnostics (CSV): crack onset, solver effort, step rate, memory ---
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
  [mem_gib]
    type = MemoryUsage
    value_type = max_process
    report_peak_value = true
    mem_units = gibibytes
  []
[]
"""


def discovery_deck_texts(case, t_fine):
    """(elasticity.i, fracture.i) of the fresh-start Discovery run."""
    elastic, fracture = cases.deck_texts(case)
    elastic = cases.substitute(elastic, restart.PRODUCTION_DT,
                               f"x = '0   {t_fine:g}'\n    y = '{restart.COARSE_DT:g} {restart.FINE_DT:g}'",
                               "discovery dt_fn")
    elastic = cases.substitute(elastic, "  end_time = 20\n", f"  end_time = {restart.END_TIME}\n",
                               "discovery end_time")
    start = elastic.index("sync_times = '")
    end = elastic.index("'", start + len("sync_times = '"))
    coarse = [f"{float(k):.1f}" for k in range(1, int(t_fine) + 1)]
    if t_fine != int(t_fine):                       # e.g. 26.5: also output where the fine steps begin
        coarse.append(f"{t_fine:.1f}")
    n_fine = int(round((restart.END_TIME - t_fine) / restart.SYNC_STEP))
    fine = [f"{t_fine + restart.SYNC_STEP * k:.1f}" for k in range(1, n_fine + 1)]
    elastic = elastic[:start] + "sync_times = '" + " ".join(coarse + fine) + elastic[end:]
    elastic = cases.substitute(elastic, "[Outputs]\n", DIAGNOSTICS + "\n[Outputs]\n  csv = true\n",
                               "discovery outputs")
    return elastic, fracture


def main():
    if len(sys.argv) != 2:
        print(__doc__)
        return 2
    tag = sys.argv[1]
    case = next(c for c in geos.CASES if geos.case_tag(c) == tag)
    t_crack, source = restart.crack_time(tag)
    t_fine = restart.fine_start(t_crack)
    case_dir = os.path.join(geos.TPB_DIR, f"code_{case['family']}", cases.case_dir_name(case) + SUFFIX)
    elastic, fracture = discovery_deck_texts(case, t_fine)
    cases.write(os.path.join(case_dir, "elasticity.i"), elastic)
    cases.write(os.path.join(case_dir, "fracture.i"), fracture)
    print(f"wrote {os.path.relpath(case_dir, geos.TPB_DIR)}/  (crack {t_crack:.1f} {source}; "
          f"dt 0.1 to {t_fine:g}, then 0.001 to {restart.END_TIME})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
