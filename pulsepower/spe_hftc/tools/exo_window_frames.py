# Extract the damage field at the end of given pulse windows (last frame at or before k*Tw) from an Exodus file into npz files.
#   python3 exo_window_frames.py OUT.e CASE_NAME OUT_DIR Tw_us k1 k2 ...
import sys, numpy as np, netCDF4
src, case, outdir, Tw = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4]); ks = [int(a) for a in sys.argv[5:]]
ds = netCDF4.Dataset(src)
t = np.asarray(ds["time_whole"][:]) * 1e6
names = [s.strip() for s in netCDF4.chartostring(ds["name_nod_var"][:])]
x, y = (ds["coordx"][:], ds["coordy"][:]) if "coordx" in ds.variables else (ds["coord"][:][0], ds["coord"][:][1])
tri = np.vstack([ds[v][:] for v in ds.variables if v.startswith("connect")]) - 1
iv = names.index("d") + 1
for k in ks:
    f = int(np.searchsorted(t, k * Tw + 1e-6, side="right") - 1)
    d = np.asarray(ds[f"vals_nod_var{iv}"][f, :], dtype=np.float32)
    np.savez_compressed(f"{outdir}/{case}_pulse{k}.npz", x=np.asarray(x), y=np.asarray(y), tri=np.asarray(tri, dtype=np.int32), t=float(t[f]) * 1e-6, d=d)
    print(case, "pulse", k, "frame", f, "t_us", round(float(t[f]), 1))
