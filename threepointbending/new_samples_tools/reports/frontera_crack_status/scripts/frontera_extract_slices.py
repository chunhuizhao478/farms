"""
frontera_extract_slices.py -- per case, the tetrahedra that cross the mid-thickness plane, for crack snapshots (read-only).

WHAT    for each new-sample case folder code_<family>/case_<family>_3d_<tag>/ on Frontera (probe and Discovery copies
        skipped), opens elasticity_exodus.e (classic netCDF, memory-mapped), takes the last complete output (the final
        record may be mid-write, so the one before it) and keeps the tetrahedra that cross z = Z_MM inside the x window
        X_WINDOW_MM.  Saves <case>.npz in OUT: xyz (mm, float32) and d (float32) of the kept nodes, tets (int32,
        renumbered), t, crack tip (highest node with d > 0.95 over the WHOLE mesh, mm), broken-node count, beam top (mm).
STEPS   1 read coordinates, connectivity (all element blocks), times, d   2 select tets   3 renumber + save
USAGE   ssh HOST 'module load python3 >/dev/null 2>&1; python3 -' < frontera_extract_slices.py
"""
import glob
import os

import numpy as np
from scipy.io import netcdf_file

P = "/scratch2/10024/zhaochun/projects/farms_cdms_04192026/threepointbending"
OUT = P + "/report_snapshots"
Z_MM = 4.0                     # mid-thickness of the 8 mm beam
X_WINDOW_MM = (6.5, 21.5)      # around the notch at x = 14 mm
D_BROKEN = 0.95


def last_frame(f):
    """Output time and d of the last complete record."""
    times = f.variables["time_whole"][:].copy()
    names = ["".join(c.decode() for c in row if c != b"").strip() for row in f.variables["name_nod_var"][:]]
    dvar = f.variables["vals_nod_var%d" % (names.index("d") + 1)]
    n = len(times) - 1
    return float(times[n - 1]), np.array(dvar[n - 1, :])


def main():
    os.makedirs(OUT, exist_ok=True)
    print("__B__")
    for path in sorted(glob.glob(P + "/code_*/case_*_3d_*/elasticity_exodus.e")):
        case = os.path.basename(os.path.dirname(path))
        if "dt0p1" in case or "discovery" in case:
            continue
        f = netcdf_file(path, "r", mmap=True)
        xyz = np.stack([np.array(f.variables[c][:]) for c in ("coordx", "coordy", "coordz")], axis=1)
        scale = 1e3 if np.abs(xyz).max() < 1.0 else 1.0
        xyz = xyz * scale
        nblk = f.dimensions["num_el_blk"]
        tets = np.concatenate([np.array(f.variables["connect%d" % (b + 1)][:]) for b in range(nblk)]) - 1
        t, d = last_frame(f)
        broken = d > D_BROKEN
        tip = float(xyz[broken, 1].max()) if broken.any() else float("nan")
        z = xyz[tets, 2]
        x = xyz[tets, 0]
        keep = (z.min(axis=1) <= Z_MM) & (z.max(axis=1) >= Z_MM) & \
               (x.max(axis=1) >= X_WINDOW_MM[0]) & (x.min(axis=1) <= X_WINDOW_MM[1])
        sel = tets[keep]
        used, inverse = np.unique(sel, return_inverse=True)
        np.savez_compressed(os.path.join(OUT, case + ".npz"), xyz=xyz[used].astype(np.float32),
                            d=d[used].astype(np.float32), tets=inverse.reshape(-1, 4).astype(np.int32),
                            t=t, tip=tip, broken=int(broken.sum()), top=float(xyz[:, 1].max()))
        print("%-62s t=%7.3f tip=%5.2f kept %6d of %7d tets" % (case, t, tip, len(sel), len(tets)))
        f.close()
    print("__E__")


if __name__ == "__main__":
    main()
