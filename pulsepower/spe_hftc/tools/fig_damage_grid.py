# Paper figure: damage field of the homogeneous-Gc cases (rows: eta 0.4 / 0.45 / 0.5) at the end of pulses 4, 8, 12, 16, 20
# (columns). Reads the npz frames written by exo_window_frames.py (last Exodus frame at or before each window end; frames are
# every 100 steps, so at most 4 us before the end, when the sample is at rest). Dashed circle through the tip of the longest
# crack; the extent is printed under each panel.
#   python3 fig_damage_grid.py OUT_STEM   -> OUT_STEM.pdf and OUT_STEM.png
import sys, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from matplotlib.patches import Circle
FR = "/private/tmp/claude-501/-Users-chunhuizhao-projects-farms-cdms/4d45f359-e932-438b-a688-5f487c854dc5/scratchpad/vista_exo/frames/"
STEM = sys.argv[1]
RH, TW = 4.15, 60.0
ROWS = [("att_homog_eta0p4_absorb", r"$\eta = 0.4$"), ("att_homog_eta0p45_absorb", r"$\eta = 0.45$"), ("att_homog_eta0p5_absorb", r"$\eta = 0.5$")]
PULSES = [4, 8, 12, 16, 20]
plt.rcParams.update({"font.size": 9, "mathtext.fontset": "dejavusans", "pdf.fonttype": 42, "ps.fonttype": 42})

fig, axes = plt.subplots(len(ROWS), len(PULSES), figsize=(10.0, 6.3), dpi=300)
pc = None
for i, (case, lab) in enumerate(ROWS):
    for j, k in enumerate(PULSES):
        ax = axes[i, j]
        z = np.load(f"{FR}{case}_pulse{k}.npz"); x, y, tri, d, t = z["x"] * 1e3, z["y"] * 1e3, z["tri"], z["d"], float(z["t"]) * 1e6
        r = np.hypot(x, y); ext = r[d >= 0.5].max() - RH if (d >= 0.5).any() else 0.0
        pc = ax.tripcolor(Triangulation(x, y, tri), d, shading="gouraud", cmap="Reds", vmin=0, vmax=1, rasterized=True)
        ax.add_patch(Circle((0, 0), RH, fill=False, ec="black", lw=0.5))
        if ext > 0: ax.add_patch(Circle((0, 0), RH + ext, fill=False, ec="black", lw=0.8, ls=(0, (4, 3)), zorder=5))
        ax.set_xlabel(f"{ext:.1f} mm", fontsize=7.5, labelpad=2)       # extent below the panel, never over the field or the circle
        ax.set_xlim(x.min(), x.max()); ax.set_ylim(y.min(), y.max()); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values(): sp.set_linewidth(0.5)
        if i == 0: ax.set_title(f"pulse {k}", fontsize=9, pad=4)
        if j == 0: ax.set_ylabel(lab, fontsize=10, labelpad=6)
        print(f"{case} pulse {k}: t = {t:.1f} us, extent {ext:.2f} mm")
fig.subplots_adjust(left=0.04, right=0.93, top=0.95, bottom=0.04, wspace=0.06, hspace=0.16)
cax = fig.add_axes([0.945, 0.15, 0.012, 0.7]); cb = fig.colorbar(pc, cax=cax); cb.set_label("$d$"); cb.ax.tick_params(labelsize=8)
fig.savefig(STEM + ".pdf"); fig.savefig(STEM + ".png", dpi=300); print("wrote", STEM + ".pdf/.png")
