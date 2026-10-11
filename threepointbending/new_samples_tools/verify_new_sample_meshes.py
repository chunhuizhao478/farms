#!/usr/bin/env python3
"""
verify_new_sample_meshes.py -- check every new-sample mesh against its STL and
against what the MOOSE decks need.

WHAT
    For each case of make_new_sample_geos.CASES whose .msh exists, prints one
    line and appends it to new_samples_tools/logs/verify_log.tsv:

      V_rel   mesh volume vs the analytic volume of the case
      d_m2s   max distance, mesh boundary nodes -> STL surface (um)
      d_s2m   max distance, STL vertices -> mesh boundary surface (um)
              (both exclude the two designed departures from the STL: the notch tip
              when NOTCH_ROUND = 1, and the throat between touching holes)
      lines   nodes on the left/right support lines and the loading line
      q_min   minimum tetrahedron mean-ratio quality (1 = regular tetrahedron)
      h_band  median edge length in the refined band above the notch (um)

STEPS
    1. read_tet_mesh()        meshio; tetrahedra, physical tags
    2. boundary_triangles()   faces used by exactly one tetrahedron
    3. stl_in_model_frame()   STL (mm) -> model frame (m); the sine-wave STLs have
                              height along z and thickness along y, so they are
                              rotated (x, y, z) -> (x, z, 8 mm - y)
    4. max_distance_um()      VTK implicit (exact point-to-triangle) distance
    5. line_node_z()          node z-coordinates on a support/loading line
    6. mean_ratio_quality()   12 (3 V)^(2/3) / sum(edge^2)

CHECKS (all must hold for PASS)
    one physical tag; no inverted tetrahedra; no orphan nodes; |V_rel| < VOL_TOL;
    d_m2s and d_s2m < the family tolerance; every BC line has >= MIN_LINE_NODES
    nodes and spans z = 0 .. 8 mm, and MOOSE's +-1e-5 m nodeset boxes catch
    exactly those nodes; q_min > Q_MIN.

USAGE
    python verify_new_sample_meshes.py [--only SUBSTRING]
"""

import argparse
import os
import sys

import meshio
import numpy as np
import pyvista as pv
import trimesh

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import make_new_sample_geos as geos  # noqa: E402

VOL_TOL = 1e-3
DIST_TOL_UM = {"2holes_asym": 10.0, "2holes_connected": 10.0, "pillars": 10.0,
               "sinewave": 40.0}   # the STL sine is a fit-point spline: <= 30 um off the sine
MIN_LINE_NODES = 9
Q_MIN = 0.05
MOOSE_BOX_HALF = 1e-5      # half width of the decks' BoundingBoxNodeSetGenerator boxes
EXACT = 1e-12              # coordinate tolerance for "on the line"
LOG_TSV = os.path.join(geos.TOOLS_DIR, "logs", "verify_log.tsv")


def read_tet_mesh(msh):
    """Points (m), tetrahedra (n, 4) and the set of physical tags of the tetrahedra."""
    m = meshio.read(msh)
    tets = np.vstack([c.data for c in m.cells if c.type == "tetra"])
    tags = set()
    for c, data in zip(m.cells, m.cell_data.get("gmsh:physical", [])):
        if c.type == "tetra":
            tags |= set(np.unique(data).tolist())
    return m.points, tets, tags


def tet_signed_volumes(p, t):
    a, b, c, d = (p[t[:, i]] for i in range(4))
    return np.einsum("ij,ij->i", b - a, np.cross(c - a, d - a)) / 6.0


def mean_ratio_quality(p, t, vol):
    """Mean-ratio shape quality in (0, 1]; 1 for a regular tetrahedron."""
    edges = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    s = sum(((p[t[:, i]] - p[t[:, j]]) ** 2).sum(1) for i, j in edges)
    return 12.0 * (3.0 * np.abs(vol)) ** (2.0 / 3.0) / s


def boundary_triangles(t):
    """Triangles (n, 3) that belong to exactly one tetrahedron."""
    faces = np.vstack([t[:, [0, 1, 2]], t[:, [0, 1, 3]], t[:, [0, 2, 3]], t[:, [1, 2, 3]]])
    key = np.sort(faces, axis=1).astype(np.int64)
    code = (key[:, 0] << 42) | (key[:, 1] << 21) | key[:, 2]
    _, first, counts = np.unique(code, return_index=True, return_counts=True)
    return faces[first[counts == 1]]


def stl_in_model_frame(case):
    """The case STL as a trimesh in the model frame (metres)."""
    m = trimesh.load_mesh(os.path.join(geos.STL_ROOT, case["stl"]), process=True)
    v = m.vertices.copy()
    if case["family"] == "sinewave":
        v = np.column_stack([v[:, 0], v[:, 2], geos.BEAM_THICKNESS_MM - v[:, 1]])
    return trimesh.Trimesh(v * 1e-3, m.faces, process=False)


def outside_designed_deviations(case, xy):
    """Mask of points outside the two regions where the model departs from the STL on purpose:
    the notch tip (round vs the STL's flat tip, when NOTCH_ROUND = 1) and, for touching holes,
    the throat that trims the zero-width material cusps (make_new_sample_geos.throat_half_height_mm)."""
    keep = np.ones(len(xy), bool)
    x0, y0 = geos.NOTCH_TIP_MM[0] * 1e-3, geos.NOTCH_TIP_MM[1] * 1e-3
    if geos.NOTCH_ROUND == 1:
        half_w = geos.NOTCH_WIDTH_MM / 2 * 1e-3
        keep &= ~((np.abs(xy[:, 0] - x0) < half_w + 1e-5)
                  & (xy[:, 1] > y0 - half_w - 1e-5) & (xy[:, 1] < y0 + 1e-5))
    if case["family"] == "2holes_connected" and geos.holes_touch(case):
        yc = (geos.NOTCH_TIP_MM[1] + case["dy"]) * 1e-3
        half_gap = (geos.CUSP_MIN_WIDTH_MM / 2 + 0.01) * 1e-3     # trimmed flat face + margin
        half_h = (geos.throat_half_height_mm(case) + 0.01) * 1e-3
        keep &= ~((np.abs(xy[:, 0] - x0) < half_gap) & (np.abs(xy[:, 1] - yc) < half_h))
    return keep


def polydata(points, triangles):
    return pv.PolyData(points, np.hstack([np.full((len(triangles), 1), 3), triangles]).ravel())


def max_distance_um(query_points, surface):
    """Largest unsigned distance from the points to a triangulated surface (um)."""
    d = pv.PolyData(query_points).compute_implicit_distance(surface)["implicit_distance"]
    return float(np.abs(d).max() * 1e6)


def line_node_z(p, x, y, half):
    """Sorted z of nodes within `half` of the line (x, y, z) in x and y."""
    on = (np.abs(p[:, 0] - x) <= half) & (np.abs(p[:, 1] - y) <= half)
    return np.sort(p[on, 2])


def verify_case(case):
    """Return (ok, row_dict) for one case."""
    msh = geos.geo_path(case)[:-4] + ".msh"
    p, t, tags = read_tet_mesh(msh)
    vol = tet_signed_volumes(p, t)
    q = mean_ratio_quality(p, t, vol)
    orphans = len(p) - len(np.unique(t))

    v_mesh = np.abs(vol).sum() * 1e9
    v_rel = (v_mesh - geos.expected_volume_mm3(case)) / geos.expected_volume_mm3(case)

    tri = boundary_triangles(t)
    bnodes = np.unique(tri)
    stl = stl_in_model_frame(case)
    keep_mesh = outside_designed_deviations(case, p[bnodes])
    keep_stl = outside_designed_deviations(case, stl.vertices)
    d_m2s = max_distance_um(p[bnodes][keep_mesh], polydata(stl.vertices, stl.faces))
    d_s2m = max_distance_um(stl.vertices[keep_stl], polydata(p, tri))

    top = geos.beam_top_mm(case) * 1e-3
    lines = {"left": (0.004, 0.0), "right": (0.024, 0.0), "load": (0.014, top)}
    line_ok, line_counts = True, []
    for name, (x, y) in lines.items():
        z_exact = line_node_z(p, x, y, EXACT)
        z_box = line_node_z(p, x, y, MOOSE_BOX_HALF)
        ok = (len(z_exact) >= MIN_LINE_NODES and len(z_box) == len(z_exact)
              and abs(z_exact[0]) < EXACT and abs(z_exact[-1] - 0.008) < EXACT)
        line_ok &= ok
        line_counts.append(len(z_exact))

    centroid = p[t].mean(axis=1)
    band = ((np.abs(centroid[:, 0] - 0.014) < 0.0008) & (centroid[:, 1] > 0.002) & (centroid[:, 1] < 0.0035))
    tb = t[band]
    edge_len = np.concatenate([np.linalg.norm(p[tb[:, i]] - p[tb[:, j]], axis=1)
                               for i, j in [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]])

    tol = DIST_TOL_UM[case["family"]]
    ok = (len(tags) == 1 and (vol < 0).sum() == 0 and orphans == 0 and abs(v_rel) < VOL_TOL
          and d_m2s < tol and d_s2m < tol and line_ok and q.min() > Q_MIN)
    row = dict(mesh=geos.mesh_basename(case), status="PASS" if ok else "FAIL", nodes=len(p), tets=len(t),
               tags=",".join(str(x) for x in sorted(tags)), inverted=int((vol < 0).sum()), orphans=orphans,
               v_rel=f"{v_rel:+.1e}", d_m2s_um=f"{d_m2s:.1f}", d_s2m_um=f"{d_s2m:.1f}", tol_um=f"{tol:g}",
               line_nodes="/".join(map(str, line_counts)), lines_ok=line_ok,
               q_min=f"{q.min():.3f}", q_lt_0p1=int((q < 0.1).sum()),
               h_band_um=f"{np.median(edge_len) * 1e6:.0f}")
    return ok, row


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", default="")
    args = ap.parse_args()
    rows, failures = [], 0
    for case in geos.CASES:
        name = geos.mesh_basename(case)
        if args.only not in name or not os.path.exists(geos.geo_path(case)[:-4] + ".msh"):
            continue
        ok, row = verify_case(case)
        failures += not ok
        rows.append(row)
        print("  ".join(f"{k}={v}" for k, v in row.items()), flush=True)
    if rows:
        os.makedirs(os.path.dirname(LOG_TSV), exist_ok=True)
        with open(LOG_TSV, "w") as fh:
            fh.write("\t".join(rows[0]) + "\n")
            for r in rows:
                fh.write("\t".join(str(v) for v in r.values()) + "\n")
    print(f"{len(rows) - failures}/{len(rows)} meshes PASS")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
