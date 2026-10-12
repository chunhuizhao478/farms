#!/usr/bin/env python3
"""
make_restart_decks.py -- turn the 16 production decks into restart decks that resolve the crack.

WHAT
    The production runs stop at end_time = 20, before any crack forms (dt probes, 2026-10-12).
    Their restart (jobs_<family>_recover/, `--recover` from the hourly checkpoint near t = 20)
    reads the same code_<family>/case_*/elasticity.i, so that file is edited in place
    (MOOSE names checkpoints and Exodus output after the input file; the production settings
    stay in git history, commit 01fb4f9).  Per case only these change:
      - dt_fn       0.1 to t = 14, 0.001 to 20 (as production), 0.1 from 20 to t_fine,
                    0.001 from t_fine on;  t_fine = crack time - FINE_LEAD, rounded down to
                    0.5 and never before 20
      - end_time    50 for every case (user, 2026-10-12: the propagation time is unknown; the
                    48 h wall time ends each run first)
      - sync_times  production list (to 20.0) continued every 0.1 to 50

    Crack times (first t with max d > 0.95 at dt = 0.1, make_dt_probe.py) were MEASURED for five
    cases and FITTED for the other eleven from their production max d at t = 15.5:
        t_crack = FIT_A + FIT_B / maxd_15p5      (5 probes, residuals <= 0.3 time units)

STEPS
    1. crack_time()       measured value, else the fit
    2. fine_start()       t_fine from the crack time
    3. restart_deck()     asserted substitutions on the production deck
    4. write the schedule to logs/restart_schedule.tsv

CHECKS
    Every substitution must match exactly once (make_new_sample_cases.substitute).

USAGE
    python make_restart_decks.py            # writes all 16 decks + the schedule
"""

import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import make_new_sample_cases as cases  # noqa: E402
import make_new_sample_geos as geos  # noqa: E402

FINE_LEAD = 1.5          # start 0.001 steps this long before the crack (hard steps begin 0.7-1.0 before)
COARSE_DT = 0.1          # accurate before the crack: probe max d within 2-3 % of production
FINE_DT = 0.001
RESTART_FROM = 20.0      # production end_time; restarts recover from the checkpoint just before it
END_TIME = 50
SYNC_STEP = 0.1
FIT_A, FIT_B = 14.40, 1.521

# production max d at t = 15.5 (Exodus scan 2026-10-12) and measured crack times (dt probes)
MAXD_15P5 = {
    "xL3.00mm_xR2.00mm_y3.6mm_radi1.0mm": 0.2011, "xL3.00mm_xR2.50mm_y3.6mm_radi1.0mm": 0.1897,
    "xL3.50mm_xR2.50mm_y3.6mm_radi1.0mm": 0.1754,
    "DL1.5mm_DR2.0mm_x1.00mm_y3.6mm_overlap": 0.1500, "DL1.5mm_DR2.0mm_x1.50mm_y3.6mm_downline": 0.1136,
    "DL1.5mm_DR2.0mm_x1.50mm_y3.6mm_upline": 0.1164, "DL1.5mm_DR2.0mm_x2.00mm_y3.6mm_downline": 0.1101,
    "DL1.5mm_DR2.0mm_x2.00mm_y3.6mm_midline": 0.1010, "DL2.0mm_DR2.0mm_x1.00mm_y3.6mm_overlap": 0.1110,
    "DL2.0mm_DR2.0mm_x1.50mm_y3.6mm_downline": 0.1329, "DL2.0mm_DR2.0mm_x2.00mm_y3.6mm_downline": 0.1097,
    "DL2.0mm_DR2.0mm_x2.00mm_y3.6mm_upline": 0.0958,
    "D1.0mm_y1.6mm_amp2.0mm_centered": 0.1760, "D1.0mm_y1.6mm_amp2.0mm_1ampoffset": 0.1800,
    "welded_mismatch_2onN": 0.1550, "welded_mismatch_3onN": 0.1528,
}
MEASURED_CRACK = {
    "xL3.00mm_xR2.00mm_y3.6mm_radi1.0mm": 21.8,          # job 7990155
    "D1.0mm_y1.6mm_amp2.0mm_centered": 22.9,             # job 7990181
    "welded_mismatch_2onN": 24.5,                        # job 7990181
    "DL1.5mm_DR2.0mm_x1.50mm_y3.6mm_downline": 28.0,     # job 7990271
    "DL2.0mm_DR2.0mm_x2.00mm_y3.6mm_upline": 30.0,       # job 7990271
}
PRODUCTION_DT = "x = '0   14'\n    y = '0.1 0.001'"


def crack_time(tag):
    """(crack time, 'measured' | 'fitted')."""
    if tag in MEASURED_CRACK:
        return MEASURED_CRACK[tag], "measured"
    return FIT_A + FIT_B / MAXD_15P5[tag], "fitted"


def fine_start(t_crack):
    """Start of the 0.001 steps: FINE_LEAD before the crack, down to a multiple of 0.5, >= 20."""
    return max(RESTART_FROM, math.floor((t_crack - FINE_LEAD) * 2) / 2)


def restart_deck(case, t_fine):
    """Production elasticity.i with the restart dt schedule, end_time and output times."""
    text, _ = cases.deck_texts(case)
    if t_fine > RESTART_FROM:
        new_dt = (f"x = '0   14    {RESTART_FROM:g}   {t_fine:g}'\n"
                  f"    y = '0.1 {FINE_DT:g} {COARSE_DT:g} {FINE_DT:g}'")
    else:                              # crack too close to 20: stay at 0.001
        new_dt = PRODUCTION_DT
    text = cases.substitute(text, PRODUCTION_DT, new_dt, "restart dt_fn") if new_dt != PRODUCTION_DT else text
    text = cases.substitute(text, "  end_time = 20\n", f"  end_time = {END_TIME}\n", "restart end_time")
    old_tail = "19.9 20.0'"
    n = int(round((END_TIME - RESTART_FROM) / SYNC_STEP))
    extra = " ".join(f"{RESTART_FROM + SYNC_STEP * k:.1f}" for k in range(1, n + 1))
    text = cases.substitute(text, old_tail, f"19.9 20.0 {extra}'", "restart sync_times")
    return text


def main():
    rows = ["case\tmaxd_15p5\tcrack_t\tsource\tfine_start\tend_time"]
    for case in geos.CASES:
        tag = geos.case_tag(case)
        t_crack, source = crack_time(tag)
        t_fine = fine_start(t_crack)
        path = os.path.join(geos.TPB_DIR, f"code_{case['family']}", cases.case_dir_name(case), "elasticity.i")
        cases.write(path, restart_deck(case, t_fine))
        rows.append(f"{tag}\t{MAXD_15P5[tag]:.4f}\t{t_crack:.1f}\t{source}\t{t_fine:.1f}\t{END_TIME}")
        print(f"{tag:45s} crack {t_crack:5.1f} ({source:8s}) -> 0.1 from 20 to {t_fine:4.1f}, "
              f"then 0.001 to {END_TIME}")
    with open(os.path.join(geos.TOOLS_DIR, "logs", "restart_schedule.tsv"), "w") as fh:
        fh.write("\n".join(rows) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
