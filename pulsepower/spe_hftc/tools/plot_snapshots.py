# 2x2 damage snapshots of the four Vista cases (latest frame of each): python3 plot_snapshots.py OUT.png
import sys, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from matplotlib.patches import Circle
S = "/private/tmp/claude-501/-Users-chunhuizhao-projects-farms-cdms/4d45f359-e932-438b-a688-5f487c854dc5/scratchpad/vista_exo/"
OUT = sys.argv[1]
INK, INK2, PAPER = "#0b0b0b", "#52514e", "#fcfcfb"; RH, RCORE = 4.15, 15.0
CASES = [("att_eta0p5_absorb", "Production: Gaussian Gc, eta 0.5, spectral split", "job 1089067 + tail 1092063"),
         ("att_homog_eta0p5_absorb", "Homogeneous Gc, eta 0.5, spectral split", "job 1090509, recovering as 1092069"),
         ("att_homog_eta0p4_absorb", "Homogeneous Gc, eta 0.4, spectral split", "job 1090509, recovering as 1092069"),
         ("att_eta0p5_voldev_absorb", "Gaussian Gc, eta 0.5, volumetric-deviatoric split", "job 1090603, stopped by the file-system event")]
def cracks_of(x, y, d):
    r = np.hypot(x, y); th = np.degrees(np.arctan2(y, x)) % 360
    m = d >= 0.5; sec = np.zeros(360)
    for a, rr in zip(th[m].astype(int), r[m]): sec[a] = max(sec[a], rr)
    ext = sec - RH; hit = ext > 2.0; cr = []; i = 0
    while i < 360:
        if hit[i]:
            j = i
            while j < 360 and hit[j]: j += 1
            cr.append((ext[i:j].max(), (i + j) / 2)); i = j
        else: i += 1
    return sorted(cr, reverse=True)
fig, axes = plt.subplots(2, 2, figsize=(11.6, 11.0), dpi=150)
fig.patch.set_facecolor(PAPER); fig.subplots_adjust(left=0.07, right=0.88, top=0.90, bottom=0.06, wspace=0.22, hspace=0.30)
for ax, (c, title, sub) in zip(axes.flat, CASES):
    z = np.load(S + c + "_last.npz"); x, y, tri, d, t = z["x"]*1e3, z["y"]*1e3, z["tri"], z["d"], float(z["t"])*1e6
    r = np.hypot(x, y); pulse = int(t // 60.0) + 1
    e5 = r[d >= 0.5].max() - RH if (d >= 0.5).any() else 0.0
    core = r < RCORE; e5c = r[(d >= 0.5) & core].max() - RH if ((d >= 0.5) & core).any() else 0.0
    cr = cracks_of(x, y, d)
    pc = ax.tripcolor(Triangulation(x, y, tri), d, shading="gouraud", cmap="Reds", vmin=0, vmax=1, rasterized=True)
    ax.add_patch(Circle((0, 0), RH, fill=False, ec=INK2, lw=0.8)); ax.add_patch(Circle((0, 0), RCORE, fill=False, ec=INK2, lw=0.8, ls=(0, (3, 3))))
    half = max(abs(x).max(), abs(y).max()); ax.set_xlim(-half, half); ax.set_ylim(-half, half); ax.set_aspect("equal")
    ax.set_facecolor(PAPER)
    m5 = d >= 0.5; far = np.argmax(np.where(m5, r, -1.0)); corner = abs(x[far]) > 24 and abs(y[far]) > 24
    # damage continuous around the hole? (fraction of 1-degree sectors with d >= 0.5 beyond 2 mm)
    th = np.degrees(np.arctan2(y, x)) % 360; sec = np.zeros(360, bool)
    for a in th[m5 & (r > RH + 2.0)].astype(int): sec[a] = True
    ring = sec.mean() > 0.9
    body = r[m5 & ~((abs(x) > 24) & (abs(y) > 24))].max() - RH if m5.any() else 0.0
    line1 = f"{t:.0f} us, in pulse {pulse}"
    if ring:
        line2 = f"damage continuous around the hole to {body:.1f} mm"
        line3 = f"(crushed zone, no discrete cracks)" + (f"\nplus a damaged spot at the lower-left corner, {e5:.1f} mm" if corner else "")
    else:
        line2 = f"farthest d >= 0.5: {body:.1f} mm from the wall" + (" (at the side face)" if body > 21.3 else "")
        line3 = f"{len(cr)} radial cracks longer than 2 mm" + (f", longest {cr[0][0]:.1f} mm" if cr else "")
    ax.text(0.02, 0.985, line1 + "\n" + line2 + "\n" + line3, transform=ax.transAxes, fontsize=8.3, color=INK, va="top",
            bbox=dict(fc=PAPER, ec="none", alpha=0.85, pad=2))
    ax.set_title(title + "\n" + sub, fontsize=9.5, color=INK, loc="left")
    ax.set_xlabel("x (mm)", color=INK2, fontsize=9); ax.set_ylabel("y (mm)", color=INK2, fontsize=9); ax.tick_params(colors=INK2, labelsize=8)
    for s in ax.spines.values(): s.set_color("#c9c8c3")
    print(f"{c:26s} t={t:7.1f} us pulse {pulse:2d}  ext(d>=0.5) {e5:5.2f} mm  core {e5c:5.2f} mm  cracks>2mm {len(cr)}: " + ", ".join(f"{a:.1f}@{b:.0f}" for a, b in cr[:8]))
cax = fig.add_axes([0.905, 0.25, 0.018, 0.5]); cb = fig.colorbar(pc, cax=cax); cb.set_label("damage d", color=INK2); cb.ax.tick_params(colors=INK2, labelsize=8)
fig.suptitle("Damage field of each Vista case at its latest frame (full 2D section, borehole 8.3 mm, dashed circle = 15 mm core)", x=0.02, ha="left", color=INK, fontsize=11.5, y=0.965)
fig.savefig(OUT, facecolor=PAPER); print("wrote", OUT)
