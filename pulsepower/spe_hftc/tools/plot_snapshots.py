# 2x3 damage snapshots of the Vista section cases (latest frame of each) plus the per-pulse table against the paper:
#   python3 plot_snapshots.py OUT.png      (reads vista_exo/<case>_last.npz and vista_csv/<case>.csv)
import sys, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from matplotlib.patches import Circle
S = "/private/tmp/claude-501/-Users-chunhuizhao-projects-farms-cdms/4d45f359-e932-438b-a688-5f487c854dc5/scratchpad/"
sys.path.insert(0, S); from pulse_table import load, windows
OUT = sys.argv[1]
INK, INK2, PAPER, EDGE = "#0b0b0b", "#52514e", "#fcfcfb", "#c9c8c3"; RH, RCORE = 4.15, 15.0
CASES = [("att_eta0p5_absorb", "Gaussian Gc, eta 0.5 (production)", "job 1089067 + tail 1092063", "Gaussian 0.5"),
         ("att_eta0p45_absorb", "Gaussian Gc, eta 0.45", "job 1092158", "Gaussian 0.45"),
         None,   # table panel
         ("att_homog_eta0p5_absorb", "Homogeneous Gc, eta 0.5", "job 1090509, recovered as 1092069", "homogeneous 0.5"),
         ("att_homog_eta0p45_absorb", "Homogeneous Gc, eta 0.45", "job 1092158", "homogeneous 0.45"),
         ("att_homog_eta0p4_absorb", "Homogeneous Gc, eta 0.4", "job 1090509, recovered as 1092069", "homogeneous 0.4")]
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
def table_panel(ax):
    ax.axis("off")
    cols = ["", "pulse 5", "pulse 11", "pulse 15", "pulse 20", "E_c (J/m)", "W_p (kJ/m)"]
    ct = [["experiment", "7.97", "13.90", "15.66", "20.45", "-", "-"]]
    for case in CASES:
        if case is None: continue
        rows, _ = load(S + "vista_csv/" + case[0] + ".csv"); w = {k: (e * 1e3, c) for k, e, c, t in windows(rows)}
        running = rows[-1]["time"] * 1e6 < 1199.0
        cells = [case[3]]
        for k in (5, 11, 15, 20):
            if k in w and w[k][1]: cells.append(f"{w[k][0]:.2f}")
            elif k in w: cells.append(f"({w[k][0]:.2f})")
            else: cells.append("-")
        ce, pw = rows[-1]["crackE"], -rows[-1]["external_work"] / 1e3
        cells += [f"({ce:.0f})" if running else f"{ce:.0f}", f"({pw:.1f})" if running else f"{pw:.1f}"]
        ct.append(cells)
    tb = ax.table(cellText=ct, colLabels=cols, loc="center", cellLoc="center", bbox=[-0.04, 0.36, 1.1, 0.46],
                  colWidths=[0.29, 0.1, 0.1, 0.1, 0.1, 0.155, 0.155])
    tb.auto_set_font_size(False); tb.set_fontsize(7.9)
    for (i, j), cell in tb.get_celld().items():
        cell.set_edgecolor(EDGE); cell.set_linewidth(0.6); cell.get_text().set_color(INK)
        if i == 0 or j == 0: cell.get_text().set_fontweight("bold"); cell.set_facecolor("#f1f0ec")
        if i == 0: cell.get_text().set_fontsize(7.3)
        else: cell.set_facecolor(PAPER)
        if j == 0 and i > 0: cell.get_text().set_ha("left"); cell.PAD = 0.04
    ax.set_title("Crack extent from the wall at the end of the pulse (mm)", fontsize=9.5, color=INK, loc="left")
    ax.text(0.0, 0.30, "Measured: Jin et al., 144 J pulses, 80 s apart.\nModel: farthest d >= 0.5 from the wall at the end of each\n"
            "60 us window; the nearest side face is 21.4 mm from the wall.\n( ) = pulse or run still in progress.\n"
            "E_c = crack energy, W_p = pressure work, per unit thickness at the\nlast frame; crack energy summed across the recover resets.\n"
            "Gaussian Gc, eta 0.4 has not been run.",
            transform=ax.transAxes, fontsize=8.2, color=INK2, va="top", linespacing=1.4)
fig, axes = plt.subplots(2, 3, figsize=(16.8, 11.2), dpi=150)
fig.patch.set_facecolor(PAPER); fig.subplots_adjust(left=0.045, right=0.915, top=0.905, bottom=0.055, wspace=0.2, hspace=0.26)
for ax, case in zip(axes.flat, CASES):
    if case is None: table_panel(ax); continue
    c, title, sub, short = case
    z = np.load(S + "vista_exo/" + c + "_last.npz"); x, y, tri, d, t = z["x"]*1e3, z["y"]*1e3, z["tri"], z["d"], float(z["t"])*1e6
    r = np.hypot(x, y); pulse = max(1, int(np.ceil(t / 60.0 - 1e-6)))   # window k covers ((k-1)*60, k*60] us
    e5 = r[d >= 0.5].max() - RH if (d >= 0.5).any() else 0.0
    core = r < RCORE; e5c = r[(d >= 0.5) & core].max() - RH if ((d >= 0.5) & core).any() else 0.0
    cr = cracks_of(x, y, d)
    pc = ax.tripcolor(Triangulation(x, y, tri), d, shading="gouraud", cmap="Reds", vmin=0, vmax=1, rasterized=True)
    ax.add_patch(Circle((0, 0), RH, fill=False, ec=INK2, lw=0.8)); ax.add_patch(Circle((0, 0), RCORE, fill=False, ec=INK2, lw=0.8, ls=(0, (3, 3))))
    half = max(abs(x).max(), abs(y).max()); ax.set_xlim(-half, half); ax.set_ylim(-half, half); ax.set_aspect("equal")
    ax.set_facecolor(PAPER)
    m5 = d >= 0.5; far = np.argmax(np.where(m5, r, -1.0)); corner = abs(x[far]) > 24 and abs(y[far]) > 24
    th = np.degrees(np.arctan2(y, x)) % 360; sec = np.zeros(360, bool)
    for a in th[m5 & (r > RH + 2.0)].astype(int): sec[a] = True
    ring = sec.mean() > 0.9
    body = r[m5 & ~((abs(x) > 24) & (abs(y) > 24))].max() - RH if m5.any() else 0.0
    line1 = f"{t:.0f} us, " + (f"end of pulse {pulse}" if abs(t - 60 * pulse) < 1.0 else f"in pulse {pulse}")
    if ring:
        line2 = f"damage continuous around the hole to {body:.1f} mm"
        line3 = f"(crushed zone, no discrete cracks)" + (f"\nplus a damaged spot at the lower-left corner, {e5:.1f} mm" if corner else "")
    else:
        line2 = f"farthest d >= 0.5: {body:.1f} mm from the wall" + (" (at the side face)" if body > 21.3 else "")
        line3 = f"{len(cr)} radial cracks longer than 2 mm" + (f", longest {cr[0][0]:.1f} mm" if cr else "")
    ax.text(0.02, 0.985, line1 + "\n" + line2 + "\n" + line3, transform=ax.transAxes, fontsize=8.3, color=INK, va="top",
            bbox=dict(fc=PAPER, ec="none", alpha=0.85, pad=2))
    ax.set_title(title + "\n" + sub + (", still running" if t < 1199.0 else ""), fontsize=9.5, color=INK, loc="left")
    ax.set_xlabel("x (mm)", color=INK2, fontsize=9); ax.set_ylabel("y (mm)", color=INK2, fontsize=9); ax.tick_params(colors=INK2, labelsize=8)
    for s in ax.spines.values(): s.set_color(EDGE)
    print(f"{c:26s} t={t:7.1f} us pulse {pulse:2d}  ext(d>=0.5) {e5:5.2f} mm  core {e5c:5.2f} mm  cracks>2mm {len(cr)}: " + ", ".join(f"{a:.1f}@{b:.0f}" for a, b in cr[:8]))
cax = fig.add_axes([0.93, 0.25, 0.013, 0.5]); cb = fig.colorbar(pc, cax=cax); cb.set_label("damage d", color=INK2); cb.ax.tick_params(colors=INK2, labelsize=8)
fig.suptitle("Damage field of each Vista section case at its latest frame (AT1, spectral split, absorbing side faces; borehole 8.3 mm, dashed circle = 15 mm core)",
             x=0.02, ha="left", color=INK, fontsize=11.5, y=0.965)
fig.savefig(OUT, facecolor=PAPER); print("wrote", OUT)
