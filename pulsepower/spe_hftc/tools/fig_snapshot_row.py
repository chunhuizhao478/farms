# Working figure: one row of damage snapshots, several cases at the end of one pulse, each with the dashed circle through the
# tip of its longest crack; prints the arm list (longest radius with d >= 0.5 per degree sector, contiguous sectors > 2 mm = one arm).
#   python3 fig_snapshot_row.py OUT_STEM PULSE case=label [case=label ...]    (reads vista_exo/frames/<case>_pulse<PULSE>.npz)
import sys, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from matplotlib.patches import Circle
FR = "/private/tmp/claude-501/-Users-chunhuizhao-projects-farms-cdms/4d45f359-e932-438b-a688-5f487c854dc5/scratchpad/vista_exo/frames/"
STEM, PULSE = sys.argv[1], int(sys.argv[2]); CASES = [a.split("=", 1) for a in sys.argv[3:]]
RH = 4.15   # borehole radius, mm
plt.rcParams.update({"font.size": 9, "mathtext.fontset": "dejavusans", "axes.linewidth": 0.8})
def arms(x, y, d):
    r = np.hypot(x, y); th = np.degrees(np.arctan2(y, x)) % 360
    m = d >= 0.5; sec = np.zeros(360)
    for a, rr in zip(th[m].astype(int), r[m]): sec[a] = max(sec[a], rr)
    ext = sec - RH; hit = ext > 2.0; out = []; i = 0
    while i < 360:
        if hit[i]:
            j = i
            while j < 360 and hit[j]: j += 1
            out.append((ext[i:j].max(), (i + j) / 2)); i = j
        else: i += 1
    return sorted(out, reverse=True)
n = len(CASES)
fig, axs = plt.subplots(1, n, figsize=(3.05 * n + 0.6, 3.4), dpi=300)
for ax, (case, label) in zip(np.atleast_1d(axs), CASES):
    z = np.load(f"{FR}{case}_pulse{PULSE}.npz"); x, y, tri, d, t = z["x"] * 1e3, z["y"] * 1e3, z["tri"], z["d"], float(z["t"]) * 1e6
    r = np.hypot(x, y); ext = r[d >= 0.5].max() - RH
    pc = ax.tripcolor(Triangulation(x, y, tri), d, shading="gouraud", cmap="Reds", vmin=0, vmax=1, rasterized=True)
    ax.add_patch(Circle((0, 0), RH, fill=False, ec="black", lw=0.6))
    ax.add_patch(Circle((0, 0), RH + ext, fill=False, ec="black", lw=0.9, ls=(0, (4, 3)), zorder=5))
    ax.set_xlim(x.min(), x.max()); ax.set_ylim(y.min(), y.max()); ax.set_aspect("equal")
    ax.set_xticks(range(-20, 21, 10)); ax.set_yticks(range(-20, 21, 10)); ax.tick_params(labelsize=8, length=3)
    ax.set_title(f"{label}\npulse {PULSE}, {ext:.1f} mm", fontsize=9, pad=4)
    A = arms(x, y, d)
    print(f"{case} pulse {PULSE} (t = {t:.1f} us): extent {ext:.2f} mm, {len(A)} arms > 2 mm: " + ", ".join(f"{e:.1f} mm @ {a:.0f} deg" for e, a in A))
cb = fig.colorbar(pc, ax=list(np.atleast_1d(axs)), fraction=0.02, pad=0.02, ticks=[0, 0.5, 1]); cb.set_label("$d$"); cb.ax.tick_params(labelsize=8)
fig.savefig(STEM + ".png", dpi=300, bbox_inches="tight"); print("wrote", STEM + ".png")
