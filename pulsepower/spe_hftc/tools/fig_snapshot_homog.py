# Paper figures: standalone damage snapshot of each homogeneous-Gc section case at its last frame (1200 us, end of pulse 20),
# with a dashed circle through the tip of the longest crack (radius = borehole radius + crack extent from the wall).
#   python3 fig_snapshot_homog.py OUT_DIR      -> OUT_DIR/snapshot_homog_eta0p4.{pdf,png} etc. (reads vista_exo/<case>_last.npz)
import sys, os, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from matplotlib.patches import Circle
from matplotlib.lines import Line2D
S = "/private/tmp/claude-501/-Users-chunhuizhao-projects-farms-cdms/4d45f359-e932-438b-a688-5f487c854dc5/scratchpad/vista_exo/"
OUT = sys.argv[1]; os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.size": 9, "mathtext.fontset": "dejavusans", "axes.linewidth": 0.8, "pdf.fonttype": 42, "ps.fonttype": 42})
RH = 4.15   # borehole radius, mm
CASES = [("att_homog_eta0p4_absorb", "eta0p4", r"$\eta = 0.4$"),
         ("att_homog_eta0p45_absorb", "eta0p45", r"$\eta = 0.45$"),
         ("att_homog_eta0p5_absorb", "eta0p5", r"$\eta = 0.5$")]
for case, tag, label in CASES:
    z = np.load(S + case + "_last.npz"); x, y, tri, d, t = z["x"] * 1e3, z["y"] * 1e3, z["tri"], z["d"], float(z["t"]) * 1e6
    r = np.hypot(x, y); ext = r[d >= 0.5].max() - RH
    fig, ax = plt.subplots(figsize=(3.6, 3.3), dpi=300)
    pc = ax.tripcolor(Triangulation(x, y, tri), d, shading="gouraud", cmap="Reds", vmin=0, vmax=1, rasterized=True)
    ax.add_patch(Circle((0, 0), RH, fill=False, ec="black", lw=0.6))
    ax.add_patch(Circle((0, 0), RH + ext, fill=False, ec="black", lw=0.9, ls=(0, (4, 3)), zorder=5))
    ax.set_xlim(x.min(), x.max()); ax.set_ylim(y.min(), y.max()); ax.set_aspect("equal")
    ax.set_xticks(range(-20, 21, 10)); ax.set_yticks(range(-20, 21, 10))
    ax.set_xlabel("x (mm)"); ax.set_ylabel("y (mm)"); ax.tick_params(labelsize=8, length=3)
    # the extent label sits above the axes, on the title line, so that nothing overlays the damage field or the circle
    ax.legend([Line2D([0], [0], color="black", lw=0.9, ls=(0, (4, 3)))], [f"{ext:.1f} mm"], loc="lower right", bbox_to_anchor=(1.0, 1.0),
              frameon=False, fontsize=8.5, handlelength=2.4, borderaxespad=0.15, borderpad=0.2, handletextpad=0.6)
    ax.set_title(label, loc="left", fontsize=10, pad=6)
    cb = fig.colorbar(pc, ax=ax, fraction=0.046, pad=0.04, ticks=[0, 0.5, 1]); cb.set_label("$d$", fontsize=9); cb.ax.tick_params(labelsize=8)
    cb.outline.set_linewidth(0.6)
    fig.tight_layout(pad=0.3)
    for ext_ in ("pdf", "png"): fig.savefig(os.path.join(OUT, f"snapshot_homog_{tag}.{ext_}"), dpi=300)
    plt.close(fig)
    print(f"{case}: t = {t:.1f} us, extent {ext:.2f} mm (circle radius {RH + ext:.2f} mm) -> snapshot_homog_{tag}.pdf/png")
