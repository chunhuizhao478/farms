#!/usr/bin/env python3
"""Triangular mesh of the 2D plane-strain section of the SPE-HFTC acrylic cube (gmsh OCC).

Geometry: rectangle [-lx/2, lx/2] x [-ly/2, ly/2] in the plane z = 0 minus a disk of
radius rb at the origin (the borehole). Physical groups (names are what MOOSE uses):
  surface : acrylic
  curves  : borehole, xmin, xmax, ymin, ymax
Element size: uniform h (the whole section is the fine zone).

Examples:
  python3 make_spe_section_2d.py --h 0.3e-3   --out spe_section_h0p3     # ~67 k triangles, l = 0.6 mm
  python3 make_spe_section_2d.py --h 0.15e-3  --out spe_section_h0p15    # ~0.27 M,       l = 0.3 mm
  python3 make_spe_section_2d.py --h 0.075e-3 --out spe_section_h0p075   # ~1.1 M,        l = 0.15 mm
  python3 make_spe_section_2d.py --h 0.0385e-3 --out spe_section_h0p0385 # ~4 M,          l = 77 um
"""
import argparse, json, math, sys, time

import gmsh
import numpy as np


def build(lx, ly, rb, h, out, msh_version="4.1", threads=8, verbose=True):
    t0 = time.time()
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1 if verbose else 0)
    gmsh.option.setNumber("General.NumThreads", threads)
    gmsh.model.add("spe_section")
    occ = gmsh.model.occ

    rect = occ.addRectangle(-lx / 2, -ly / 2, 0.0, lx, ly)
    disk = occ.addDisk(0.0, 0.0, 0.0, rb, rb)
    occ.cut([(2, rect)], [(2, disk)], removeObject=True, removeTool=True)
    occ.synchronize()

    surfs = gmsh.model.getEntities(2)
    if len(surfs) != 1:
        raise RuntimeError(f"expected 1 surface after cut, got {len(surfs)}")
    stag = surfs[0][1]
    curves = gmsh.model.getEntities(1)
    if len(curves) != 5:
        raise RuntimeError(f"expected 5 curves after cut, got {len(curves)}")

    tol = 1e-6
    names = {}
    for (_, ctag) in curves:
        cx, cy, _ = occ.getCenterOfMass(1, ctag)
        if abs(cx + lx / 2) < tol:
            names[ctag] = "xmin"
        elif abs(cx - lx / 2) < tol:
            names[ctag] = "xmax"
        elif abs(cy + ly / 2) < tol:
            names[ctag] = "ymin"
        elif abs(cy - ly / 2) < tol:
            names[ctag] = "ymax"
        elif abs(cx) < tol and abs(cy) < tol:
            names[ctag] = "borehole"
        else:
            raise RuntimeError(f"unclassified curve {ctag}: com=({cx},{cy})")
    expected = {"borehole", "xmin", "xmax", "ymin", "ymax"}
    if set(names.values()) != expected:
        raise RuntimeError(f"curve names {sorted(names.values())} != {sorted(expected)}")

    phys = {}
    for ctag, nm in names.items():
        phys[nm] = gmsh.model.addPhysicalGroup(1, [ctag], name=nm)
    phys["acrylic"] = gmsh.model.addPhysicalGroup(2, [stag], name="acrylic")

    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.MeshSizeMin", h)
    gmsh.option.setNumber("Mesh.MeshSizeMax", h)
    gmsh.option.setNumber("Mesh.Algorithm", 6)      # Frontal-Delaunay
    gmsh.option.setNumber("Mesh.Optimize", 1)
    gmsh.option.setNumber("Mesh.ElementOrder", 1)
    gmsh.option.setNumber("Mesh.MshFileVersion", float(msh_version))
    gmsh.option.setNumber("Mesh.Binary", 0)
    gmsh.model.mesh.generate(2)
    t_mesh = time.time() - t0

    node_tags, _, _ = gmsh.model.mesh.getNodes()
    tri_tags, _ = gmsh.model.mesh.getElementsByType(2)
    per_curve, lengths = {}, {}
    for nm in expected:
        ents = gmsh.model.getEntitiesForPhysicalGroup(1, phys[nm])
        n, L = 0, 0.0
        for e in ents:
            _, etags, _ = gmsh.model.mesh.getElements(1, e)
            n += sum(len(x) for x in etags)
            L += occ.getMass(1, e)
        per_curve[nm] = n
        lengths[nm] = L
    area = occ.getMass(2, stag)
    q = np.array(gmsh.model.mesh.getElementQualities(tri_tags, "gamma"))

    msh_path = out + ".msh"
    gmsh.write(msh_path)
    summary = {
        "file": msh_path,
        "args": dict(lx=lx, ly=ly, rb=rb, h=h, msh_version=msh_version, threads=threads),
        "physical_ids": phys,
        "n_nodes": int(len(node_tags)),
        "n_tris": int(len(tri_tags)),
        "lines_per_curve": per_curve,
        "lengths_m": lengths,
        "length_borehole_expected_m": 2 * math.pi * rb,
        "area_m2": area,
        "area_expected_m2": lx * ly - math.pi * rb ** 2,
        "gamma_min": float(q.min()),
        "gamma_mean": float(q.mean()),
        "gamma_below_0p1": int((q < 0.1).sum()),
        "mesh_seconds": round(t_mesh, 1),
        "total_seconds": round(time.time() - t0, 1),
        "gmsh_version": gmsh.__version__,
    }
    with open(out + ".json", "w") as f:
        json.dump(summary, f, indent=2)
    gmsh.finalize()
    return summary


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--lx", type=float, default=51.11e-3)
    p.add_argument("--ly", type=float, default=51.92e-3)
    p.add_argument("--rb", type=float, default=4.15e-3, help="borehole radius (m); 8.3 mm diameter")
    p.add_argument("--h", type=float, required=True, help="uniform element size (m)")
    p.add_argument("--out", required=True, help="output basename (writes <out>.msh and <out>.json)")
    p.add_argument("--msh-version", default="4.1")
    p.add_argument("--threads", type=int, default=8)
    a = p.parse_args(argv)
    s = build(a.lx, a.ly, a.rb, a.h, a.out, msh_version=a.msh_version, threads=a.threads)
    print(json.dumps(s, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
