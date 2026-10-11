#!/usr/bin/env python3
"""
mesh_new_samples.py -- mesh the new-sample .geo files with gmsh, one at a time,
under a memory kill-switch (36 GB laptop).

WHAT
    For every case of make_new_sample_geos.CASES (or those matching --only),
    runs `gmsh -3 <geo> -o <msh>` (MSH 2.2 is set inside each .geo) and appends
    one row per case to new_samples_tools/logs/mesh_log.tsv: wall time, peak
    RSS, lowest system free memory, nodes, tetrahedra, file size.

STEPS
    1. run_guarded()       start a command; every SAMPLE_S seconds read its RSS (ps)
                           and the system-wide free memory (memory_pressure -Q);
                           SIGKILL it when RSS > RSS_LIMIT_GB or free < FREE_MIN_PCT
       run_gmsh_guarded()  gmsh -3 under run_guarded(); a failed run leaves no .msh
    2. msh22_counts()      node and tetrahedron counts from the written file

CHECKS
    A case fails when gmsh exits non-zero, is killed by the guard, or the .msh
    holds no tetrahedra.  Existing meshes are skipped unless --force.

USAGE
    python mesh_new_samples.py [--only SUBSTRING] [--force]
"""

import argparse
import os
import re
import signal
import subprocess
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import make_new_sample_geos as geos  # noqa: E402

GMSH_EXE = "/Users/chunhuizhao/miniforge/envs/pythonenv/bin/gmsh"
RSS_LIMIT_GB = 20.0      # kill gmsh above this resident size
FREE_MIN_PCT = 12.0      # kill gmsh when system-wide free memory drops below this
SAMPLE_S = 2.0
LOG_DIR = os.path.join(geos.TOOLS_DIR, "logs")
SUMMARY_TSV = os.path.join(LOG_DIR, "mesh_log.tsv")
TET_TYPE = 4             # gmsh element type of a 4-node tetrahedron


def process_rss_gb(pid):
    """Resident set size of one process in GB (0 when it is gone)."""
    out = subprocess.run(["ps", "-o", "rss=", "-p", str(pid)], capture_output=True, text=True).stdout.strip()
    return int(out) / 1024 ** 2 if out else 0.0


def system_free_pct():
    """System-wide memory free percentage as reported by macOS memory_pressure."""
    out = subprocess.run(["memory_pressure", "-Q"], capture_output=True, text=True).stdout
    m = re.search(r"free percentage:\s*(\d+)%", out)
    return float(m.group(1)) if m else 100.0


def run_guarded(cmd, log_path, cwd, env=None):
    """Run `cmd` under the memory guard; return (status, wall_s, peak_rss_gb, min_free_pct).

    status is "ok", "KILLED (...)" when the guard fired, or "exit N" for a non-zero exit.
    """
    t0 = time.time()
    peak_rss, min_free = 0.0, 100.0
    with open(log_path, "w") as log:
        proc = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT, cwd=cwd, env=env)
        status = "ok"
        while proc.poll() is None:
            rss, free = process_rss_gb(proc.pid), system_free_pct()
            peak_rss, min_free = max(peak_rss, rss), min(min_free, free)
            if rss > RSS_LIMIT_GB or free < FREE_MIN_PCT:
                proc.send_signal(signal.SIGKILL)
                proc.wait()
                status = f"KILLED (rss {rss:.1f} GB, free {free:.0f} %)"
                break
            time.sleep(SAMPLE_S)
    if status == "ok" and proc.returncode != 0:
        status = f"exit {proc.returncode}"
    return status, time.time() - t0, peak_rss, min_free


def run_gmsh_guarded(geo, msh, log_path):
    """Mesh one .geo under the guard; a failed or killed run leaves no partial .msh."""
    result = run_guarded([GMSH_EXE, "-3", geo, "-o", msh], log_path, cwd=os.path.dirname(geo))
    if result[0] != "ok" and os.path.exists(msh):
        os.remove(msh)
    return result


def msh22_counts(msh):
    """(nodes, tetrahedra, other elements) of an ASCII MSH 2.2 file."""
    n_nodes = n_tets = n_other = 0
    with open(msh) as fh:
        section = None
        for line in fh:
            if line.startswith("$"):
                section = line.strip()
                header = True
                continue
            if section == "$Nodes" and header:
                n_nodes, header = int(line), False
            elif section == "$Elements":
                if header:
                    header = False
                    continue
                if int(line.split(maxsplit=2)[1]) == TET_TYPE:
                    n_tets += 1
                else:
                    n_other += 1
    return n_nodes, n_tets, n_other


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", default="", help="mesh only cases whose name contains this substring")
    ap.add_argument("--force", action="store_true", help="re-mesh even if the .msh exists")
    args = ap.parse_args()

    os.makedirs(LOG_DIR, exist_ok=True)
    new_file = not os.path.exists(SUMMARY_TSV)
    failures = 0
    with open(SUMMARY_TSV, "a") as summary:
        if new_file:
            summary.write("mesh\tstatus\twall_s\tpeak_rss_gb\tmin_free_pct\tnodes\ttets\tother_elems\tmsh_mb\n")
        for case in geos.CASES:
            name = geos.mesh_basename(case)
            if args.only not in name:
                continue
            geo = geos.geo_path(case)
            msh = geo[:-4] + ".msh"
            if os.path.exists(msh) and not args.force:
                print(f"skip {name} (exists)", flush=True)
                continue
            print(f"meshing {name} ...", flush=True)
            status, wall, rss, free = run_gmsh_guarded(geo, msh, os.path.join(LOG_DIR, name + ".gmsh.log"))
            nodes = tets = other = 0
            size_mb = 0.0
            if status == "ok":
                nodes, tets, other = msh22_counts(msh)
                size_mb = os.path.getsize(msh) / 1e6
                if tets == 0:
                    status = "no tetrahedra"
            failures += status != "ok"
            row = f"{name}\t{status}\t{wall:.0f}\t{rss:.2f}\t{free:.0f}\t{nodes}\t{tets}\t{other}\t{size_mb:.0f}"
            summary.write(row + "\n")
            summary.flush()
            print(f"  {status}: {wall:.0f} s, peak RSS {rss:.2f} GB, min free {free:.0f} %, "
                  f"{nodes} nodes, {tets} tets, {size_mb:.0f} MB", flush=True)
    print("DONE" if failures == 0 else f"DONE with {failures} failures", flush=True)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
