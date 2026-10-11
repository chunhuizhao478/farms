#!/usr/bin/env python3
"""
export_vtk_for_paraview.py -- convert every new-sample .msh to a legacy .vtk for ParaView.

WHAT
    ParaView does not read gmsh .msh files out of the box; meshfile_angles/ keeps a
    .vtk beside each .msh for that reason.  This writes meshfile_<family>/<mesh>.vtk
    (binary legacy VTK, tetrahedra only) with `gmsh <msh> -save -format vtk -bin`,
    one file at a time under the memory guard of mesh_new_samples.py.

CHECKS
    gmsh exits 0 and the .vtk exists; otherwise the partial file is removed.

USAGE
    python export_vtk_for_paraview.py [--only SUBSTRING]
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import make_new_sample_geos as geos  # noqa: E402
import mesh_new_samples as mesher  # noqa: E402


def export_one(case):
    """Write <mesh>.vtk next to <mesh>.msh; return (ok, status, size_mb)."""
    msh = geos.geo_path(case)[:-4] + ".msh"
    vtk = msh[:-4] + ".vtk"
    log = os.path.join(mesher.LOG_DIR, geos.mesh_basename(case) + ".vtk.log")
    cmd = [mesher.GMSH_EXE, msh, "-save", "-format", "vtk", "-bin", "-o", vtk]
    status, _, _, _ = mesher.run_guarded(cmd, log, cwd=os.path.dirname(msh))
    ok = status == "ok" and os.path.exists(vtk)
    if not ok and os.path.exists(vtk):
        os.remove(vtk)
    return ok, status, (os.path.getsize(vtk) / 1e6 if ok else 0.0)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", default="")
    args = ap.parse_args()
    failures = 0
    for case in geos.CASES:
        if args.only not in geos.mesh_basename(case):
            continue
        ok, status, size_mb = export_one(case)
        failures += not ok
        print(f"{'OK ' if ok else 'BAD'} {geos.mesh_basename(case)}.vtk  {status}  {size_mb:.0f} MB", flush=True)
    print("DONE" if failures == 0 else f"DONE with {failures} failures", flush=True)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
