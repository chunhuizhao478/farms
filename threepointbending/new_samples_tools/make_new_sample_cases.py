#!/usr/bin/env python3
"""
make_new_sample_cases.py -- MOOSE decks and Frontera jobs for the new samples,
copied from code_2holes/ and jobs_2holes/ (jobs_angles_recover/ for restarts).

WHAT
    For each case of make_new_sample_geos.CASES writes
        code_<family>/case_<family>_3d_<tag>/elasticity.i, fracture.i
        jobs_<family>/job_<tag>.sbatch            + submit_all.sh
        jobs_<family>_recover/job_<tag>_recover.sbatch + submit_all.sh
    The decks are the code_2holes decks with ONLY the mesh path changed; the
    pillar family (top face at y = 8.5 mm) also moves the loading line, its
    centre node and the elastic pad under the load up by 0.5 mm.  Physics,
    solver and output settings are untouched.

STEPS
    1. substitute()        replace text that must occur exactly once (asserted)
    2. deck_texts()        elasticity.i / fracture.i of one case
    3. job_text()          sbatch from jobs_2holes (10 N, 400 ranks, 48 h, split mesh)
    4. recover_job_text()  sbatch from jobs_angles_recover (--recover)
    5. submit_all_text()   loop over the folder's job files

CHECKS
    Every substitution is asserted to match exactly once in its template, so a
    template change cannot silently leave a stale mesh path or loading height.

USAGE
    python make_new_sample_cases.py
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import make_new_sample_geos as geos  # noqa: E402

TPB = geos.TPB_DIR
DECK_TEMPLATE_DIR = os.path.join(TPB, "code_2holes", "case_2holes_3d_x1.5mm_y1.6mm_radi0.5mm")
TEMPLATE_MESH = "'../../meshfile_2holes/mesh_2holes_3d_x1.50mm_y1.6mm_radi0.5mm.msh'"
JOB_TEMPLATE = os.path.join(TPB, "jobs_2holes", "job_x1.5mm_y1.6mm_radi0.5mm.sbatch")
TEMPLATE_JOB_NAME = "3dptrbend_2holes_x1.5mm_y1.6mm_radi0.5mm"
RECOVER_TEMPLATE = os.path.join(TPB, "jobs_angles_recover", "job_x1.50mm_y1.6mm_recover.sbatch")
TEMPLATE_RECOVER_NAME = "3dptrbend_whole_x1.50mm_y1.6mm_recover"

# pillar family: loading line on the top face at y = 8.5 mm instead of 8.0 mm
PILLAR_TOP_SUBSTITUTIONS = [
    ("bottom_left = '0.0139 0.007 0'\n    top_right = '0.0141 0.008 1'",
     "bottom_left = '0.0139 0.0075 0'\n    top_right = '0.0141 0.0085 1'"),
    ("# Top loading LINE: x = 0.014, y = 0.008, z = 0 to 0.008",
     "# Top loading LINE: x = 0.014, y = 0.0085 (top of the upper block), z = 0 to 0.008"),
    ("bottom_left = '${fparse loading_x - 1e-5} 0.0079 -1e-5'\n    top_right = '${fparse loading_x + 1e-5} 0.0081 0.009'",
     "bottom_left = '${fparse loading_x - 1e-5} 0.0084 -1e-5'\n    top_right = '${fparse loading_x + 1e-5} 0.0086 0.009'"),
    ("coord = '${loading_x} 0.008 ${z_center}'",
     "coord = '${loading_x} 0.0085 ${z_center}'"),
]


def read(path):
    with open(path) as fh:
        return fh.read()


def write(path, text, executable=False):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)
    if executable:
        os.chmod(path, 0o755)


def substitute(text, old, new, where):
    """Replace `old` by `new`; `old` must occur exactly once."""
    n = text.count(old)
    assert n == 1, f"{where}: expected exactly one {old!r}, found {n}"
    return text.replace(old, new)


def case_dir_name(case):
    return f"case_{case['family']}_3d_{geos.case_tag(case)}"


def deck_texts(case):
    """(elasticity.i, fracture.i) for one case."""
    mesh = f"'../../meshfile_{case['family']}/{geos.mesh_basename(case)}.msh'"
    elastic = substitute(read(os.path.join(DECK_TEMPLATE_DIR, "elasticity.i")), TEMPLATE_MESH, mesh, "elasticity.i")
    fracture = substitute(read(os.path.join(DECK_TEMPLATE_DIR, "fracture.i")), TEMPLATE_MESH, mesh, "fracture.i")
    if case["family"] == "pillars":
        for old, new in PILLAR_TOP_SUBSTITUTIONS:
            elastic = substitute(elastic, old, new, "elasticity.i (pillars)")
    return elastic, fracture


def job_text(case):
    tag, fam = geos.case_tag(case), case["family"]
    t = read(JOB_TEMPLATE)
    assert t.count(TEMPLATE_JOB_NAME) == 3, "job template: expected the job name 3 times (-J, -o, -e)"
    t = t.replace(TEMPLATE_JOB_NAME, f"3dptrbend_{fam}_{tag}")
    t = substitute(t, 'CASE_NAME="x1.5mm_y1.6mm_radi0.5mm"', f'CASE_NAME="{tag}"', "job")
    t = substitute(t, "code_2holes/case_2holes_3d_${CASE_NAME}", f"code_{fam}/case_{fam}_3d_${{CASE_NAME}}", "job")
    return t


def recover_job_text(case):
    tag, fam = geos.case_tag(case), case["family"]
    t = read(RECOVER_TEMPLATE)
    assert t.count(TEMPLATE_RECOVER_NAME) == 3, "recover template: expected the job name 3 times"
    t = t.replace(TEMPLATE_RECOVER_NAME, f"3dptrbend_{fam}_{tag}_recover")
    t = substitute(t, 'CASE_NAME="x1.50mm_y1.6mm"', f'CASE_NAME="{tag}"', "recover job")
    t = substitute(t, "code_angles/case_whole_3d_${CASE_NAME}", f"code_{fam}/case_{fam}_3d_${{CASE_NAME}}", "recover job")
    return t


def submit_all_text(family, recover):
    what = "recover jobs" if recover else "jobs"
    pattern = "job_*_recover.sbatch" if recover else "job_*.sbatch"
    return f"""#!/bin/bash
# Submit all {geos.FAMILY_TITLES[family]} {what}.
# Jobs run in this folder: their ../../farms-opt and output files are relative to it.

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR" || exit 1

for job in "$SCRIPT_DIR"/{pattern}; do
  echo "Submitting $(basename "$job")"
  sbatch "$job"
done
"""


def main():
    families = []
    for case in geos.CASES:
        fam, tag = case["family"], geos.case_tag(case)
        elastic, fracture = deck_texts(case)
        case_dir = os.path.join(TPB, f"code_{fam}", case_dir_name(case))
        write(os.path.join(case_dir, "elasticity.i"), elastic)
        write(os.path.join(case_dir, "fracture.i"), fracture)
        write(os.path.join(TPB, f"jobs_{fam}", f"job_{tag}.sbatch"), job_text(case))
        write(os.path.join(TPB, f"jobs_{fam}_recover", f"job_{tag}_recover.sbatch"), recover_job_text(case))
        print(f"wrote code_{fam}/{case_dir_name(case)}/, jobs_{fam}/job_{tag}.sbatch (+ recover)")
        if fam not in families:
            families.append(fam)
    for fam in families:
        write(os.path.join(TPB, f"jobs_{fam}", "submit_all.sh"), submit_all_text(fam, False), executable=True)
        write(os.path.join(TPB, f"jobs_{fam}_recover", "submit_all.sh"), submit_all_text(fam, True), executable=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
