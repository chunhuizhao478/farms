#!/usr/bin/env python3
"""Tetrahedral mesh of the SPE-HFTC acrylic cube with a through borehole (gmsh OCC).

Geometry: box [-lx/2, lx/2] x [-ly/2, ly/2] x [0, lz] minus a cylinder of radius rb
along z through (0, 0). Physical groups (names are what MOOSE uses):
  volume   : acrylic
  surfaces : borehole, xmin, xmax, ymin, ymax, zmin, zmax
Element size: h_fine inside the cylinder r <= r_fine (full height), h_coarse
outside, smooth tanh transition of width `transition`.

Example (production, ~17 M tets):
  python3 make_spe_cube_mesh.py --h-fine 0.3e-3 --h-coarse 1.0e-3 --r-fine 26e-3 --out spe_cube_h0p3
Example (smoke, ~0.17 M tets):
  python3 make_spe_cube_mesh.py --h-fine 1.5e-3 --h-coarse 1.5e-3 --out spe_cube_smoke
"""
import argparse, json, math, sys, time

import gmsh
import numpy as np


def build(lx, ly, lz, rb, h_fine, h_coarse, r_fine, out, msh_version="4.1",
          threads=8, transition=3e-3, optimize=True, verbose=True):
    t0 = time.time()
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1 if verbose else 0)
    gmsh.option.setNumber("General.NumThreads", threads)
    gmsh.model.add("spe_cube")
    occ = gmsh.model.occ

    # ---- geometry ----------------------------------------------------------
    box = occ.addBox(-lx / 2, -ly / 2, 0.0, lx, ly, lz)
    cyl = occ.addCylinder(0.0, 0.0, -lz, 0.0, 0.0, 3 * lz, rb)  # longer than the box
    occ.cut([(3, box)], [(3, cyl)], removeObject=True, removeTool=True)
    occ.synchronize()

    vols = gmsh.model.getEntities(3)
    if len(vols) != 1:
        raise RuntimeError(f"expected 1 volume after cut, got {len(vols)}")
    vtag = vols[0][1]
    surfs = gmsh.model.getEntities(2)
    if len(surfs) != 7:
        centres = [occ.getCenterOfMass(2, s[1]) for s in surfs]
        raise RuntimeError(f"expected 7 surfaces after cut, got {len(surfs)}: {centres}")

    # ---- classify faces by geometry, never by OCC tag ----------------------
    tol = 1e-6
    names = {}
    for (_, stag) in surfs:
        cx, cy, cz = occ.getCenterOfMass(2, stag)
        x0, y0, z0, x1, y1, z1 = gmsh.model.getBoundingBox(2, stag)
        if abs(cx + lx / 2) < tol:
            names[stag] = "xmin"
        elif abs(cx - lx / 2) < tol:
            names[stag] = "xmax"
        elif abs(cy + ly / 2) < tol:
            names[stag] = "ymin"
        elif abs(cy - ly / 2) < tol:
            names[stag] = "ymax"
        elif abs(cz) < tol and (x1 - x0) > 3 * rb:
            names[stag] = "zmin"
        elif abs(cz - lz) < tol and (x1 - x0) > 3 * rb:
            names[stag] = "zmax"
        elif abs((x1 - x0) - 2 * rb) < 1e-6 and abs((y1 - y0) - 2 * rb) < 1e-6:
            names[stag] = "borehole"
        else:
            raise RuntimeError(f"unclassified surface {stag}: com=({cx},{cy},{cz}) bbox x-extent={x1-x0}")
    expected = {"borehole", "xmin", "xmax", "ymin", "ymax", "zmin", "zmax"}
    if set(names.values()) != expected:
        raise RuntimeError(f"surface names {sorted(names.values())} != {sorted(expected)}")

    phys = {}
    for stag, nm in names.items():
        phys[nm] = gmsh.model.addPhysicalGroup(2, [stag], name=nm)
    phys["acrylic"] = gmsh.model.addPhysicalGroup(3, [vtag], name="acrylic")

    # ---- size field --------------------------------------------------------
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.MeshSizeMin", min(h_fine, h_coarse))
    gmsh.option.setNumber("Mesh.MeshSizeMax", max(h_fine, h_coarse))
    if abs(h_coarse - h_fine) > 1e-12:
        w = transition / 2.0
        fld = gmsh.model.mesh.field.add("MathEval")
        expr = (f"{h_fine:.6e} + ({h_coarse:.6e} - {h_fine:.6e})"
                f"*(0.5 + 0.5*tanh((sqrt(x*x + y*y) - {r_fine:.6e})/{w:.6e}))")
        gmsh.model.mesh.field.setString(fld, "F", expr)
        gmsh.model.mesh.field.setAsBackgroundMesh(fld)
    else:
        gmsh.option.setNumber("Mesh.MeshSizeMin", h_fine)
        gmsh.option.setNumber("Mesh.MeshSizeMax", h_fine)

    # ---- meshing -----------------------------------------------------------
    gmsh.option.setNumber("Mesh.Algorithm", 6)      # Frontal-Delaunay (surfaces)
    gmsh.option.setNumber("Mesh.Algorithm3D", 10)   # HXT (parallel)
    gmsh.option.setNumber("Mesh.Optimize", 1 if optimize else 0)
    gmsh.option.setNumber("Mesh.OptimizeNetgen", 0)
    gmsh.option.setNumber("Mesh.ElementOrder", 1)
    gmsh.option.setNumber("Mesh.MshFileVersion", float(msh_version))
    gmsh.option.setNumber("Mesh.Binary", 0)
    gmsh.model.mesh.generate(3)
    t_mesh = time.time() - t0

    # ---- statistics --------------------------------------------------------
    node_tags, _, _ = gmsh.model.mesh.getNodes()
    tet_tags, _ = gmsh.model.mesh.getElementsByType(4)
    tri_tags, _ = gmsh.model.mesh.getElementsByType(2)
    per_surface = {}
    areas = {}
    for nm in expected:
        ents = gmsh.model.getEntitiesForPhysicalGroup(2, phys[nm])
        n = 0
        a = 0.0
        for e in ents:
            _, etags, _ = gmsh.model.mesh.getElements(2, e)
            n += sum(len(x) for x in etags)
            a += occ.getMass(2, e)
        per_surface[nm] = n
        areas[nm] = a
    volume = occ.getMass(3, vtag)
    q = np.array(gmsh.model.mesh.getElementQualities(tet_tags, "gamma"))

    msh_path = out + ".msh"
    gmsh.write(msh_path)
    t_total = time.time() - t0

    summary = {
        "file": msh_path,
        "args": dict(lx=lx, ly=ly, lz=lz, rb=rb, h_fine=h_fine, h_coarse=h_coarse,
                     r_fine=r_fine, transition=transition, msh_version=msh_version,
                     threads=threads),
        "physical_ids": phys,
        "n_nodes": int(len(node_tags)),
        "n_tets": int(len(tet_tags)),
        "n_tris": int(len(tri_tags)),
        "tris_per_surface": per_surface,
        "areas_m2": areas,
        "area_borehole_expected_m2": 2 * math.pi * rb * lz,
        "volume_m3": volume,
        "volume_expected_m3": lx * ly * lz - math.pi * rb ** 2 * lz,
        "gamma_min": float(q.min()),
        "gamma_mean": float(q.mean()),
        "gamma_below_0p1": int((q < 0.1).sum()),
        "mesh_seconds": round(t_mesh, 1),
        "total_seconds": round(t_total, 1),
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
    p.add_argument("--lz", type=float, default=49.14e-3)
    p.add_argument("--rb", type=float, default=4.15e-3, help="borehole radius (m); 8.3 mm diameter")
    p.add_argument("--h-fine", type=float, required=True, help="element size inside r<=r_fine (m)")
    p.add_argument("--h-coarse", type=float, required=True, help="element size outside (m)")
    p.add_argument("--r-fine", type=float, default=26e-3, help="radius of the fine zone (m)")
    p.add_argument("--transition", type=float, default=3e-3, help="width of the fine/coarse transition (m)")
    p.add_argument("--out", required=True, help="output basename (writes <out>.msh and <out>.json)")
    p.add_argument("--msh-version", default="4.1")
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--no-optimize", action="store_true")
    a = p.parse_args(argv)
    s = build(a.lx, a.ly, a.lz, a.rb, a.h_fine, a.h_coarse, a.r_fine, a.out,
              msh_version=a.msh_version, threads=a.threads, transition=a.transition,
              optimize=not a.no_optimize)
    print(json.dumps(s, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
