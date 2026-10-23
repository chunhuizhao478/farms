# Crack extent at the end of each pulse for every Vista section case vs the acrylic-cube measurements:
#   python3 plot_extent_vs_pulse_allcases.py OUT.png      (reads vista_csv/*.csv; colour = eta, line style/marker = Gc model)
import sys, datetime
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
S = "/private/tmp/claude-501/-Users-chunhuizhao-projects-farms-cdms/4d45f359-e932-438b-a688-5f487c854dc5/scratchpad/"
sys.path.insert(0, S); from pulse_table import load, windows
OUT = sys.argv[1]
INK, INK2, GRID, PAPER = "#0b0b0b", "#52514e", "#e6e5e1", "#fcfcfb"
BLUE, MAG, GOLD, GREEN = "#2a78d6", "#c4379a", "#b8860b", "#15995f"   # validated palette (dataviz validator, light mode)
def series(path, key, partial=True, nmax=20):
    rows, _ = load(path); w = windows(rows, key)[:nmax]
    return [(k, e * 1e3) for k, e, c, t in w if c], ([(k, e * 1e3, t) for k, e, c, t in w if not c] if partial else [])
CASES = [  # file, key, label, colour, linestyle, marker, show partial, nmax
    ("att_eta0p5_absorb.csv", "fracture_extent", "Gaussian Gc, eta 0.5 (production)", BLUE, "-", "o", True, 20),
    ("att_eta0p45_absorb.csv", "fracture_extent", "Gaussian Gc, eta 0.45", MAG, "-", "o", True, 20),
    ("att_homog_eta0p5_absorb.csv", "fracture_extent", "homogeneous Gc, eta 0.5", BLUE, "--", "^", True, 20),
    ("att_homog_eta0p45_absorb.csv", "fracture_extent", "homogeneous Gc, eta 0.45", MAG, "--", "^", True, 20),
    ("att_homog_eta0p4_absorb.csv", "fracture_extent", "homogeneous Gc, eta 0.4", GOLD, "--", "^", True, 20),
    ("att_eta0p5_free.csv", "fracture_extent_core", "Gaussian Gc, eta 0.5, free side faces (core measure, caps at 10.85 mm; 3 pulses)", GREEN, ":", "D", False, 3),
]
paper = [(5, 7.97), (11, 13.90), (15, 15.66), (20, 20.45)]
fig, ax = plt.subplots(figsize=(9.4, 6.0), dpi=150)
fig.patch.set_facecolor(PAPER); ax.set_facecolor(PAPER)
ax.axhline(21.4, color=INK2, lw=0.9, ls=(0, (2, 3)), zorder=1)
ax.text(0.4, 21.4 + 0.3, "nearest side face: 21.4 mm from the wall", color=INK2, fontsize=8.5, va="bottom", ha="left")
handles, labels, ends = [], [], {}
for f, key, lab, col, ls, mk, show_part, nmax in CASES:
    done, part = series(S + "vista_csv/" + f, key, show_part, nmax)
    ms = 6.5 if mk != "^" else 7
    h, = ax.plot([p for p, _ in done], [e for _, e in done], color=col, lw=2, ls=ls, marker=mk, ms=ms, mec=PAPER, mew=1, zorder=4)
    if part:
        k, e, t = part[0]
        ax.plot([done[-1][0], k], [done[-1][1], e], color=col, lw=1.2, ls=(0, (1, 2)), zorder=3)
        ax.plot([k], [e], ls="none", marker=mk, ms=ms, mfc=PAPER, mec=col, mew=1.6, zorder=6)
        lab += f" (pulse {k} in progress, {t:.0f} us)"
    handles.append(h); labels.append(lab); ends[f] = (done, part)
    print(f"{lab}: " + ", ".join(f"{p}: {e:.2f}" for p, e in done) + (f"; partial {part[0][0]}: {part[0][1]:.2f}" if part else ""))
hp, = ax.plot([p for p, _ in paper], [e for _, e in paper], color=INK, lw=1.5, marker="s", ms=7, mec=PAPER, mew=1, zorder=5)
for p, e in paper: ax.text(p, e - 0.9, f"{e:.1f}", color=INK2, fontsize=8, ha="center", va="top")
handles.append(hp); labels.append("experiment (Jin et al.): 144 J pulses, 80 s apart")
# direct labels (text in ink tokens; the coloured mark beside the curve carries identity)
ax.text(11.2, 22.5, "Gaussian 0.5 (production)", color=BLUE, fontsize=8.8, ha="left", va="bottom", fontweight="bold")
ax.text(20.35, 23.3, "homogeneous 0.5", color=BLUE, fontsize=8.8, ha="left", va="center", fontweight="bold")
ax.text(12.6, 9.6, "homogeneous 0.4", color=GOLD, fontsize=8.8, ha="left", va="top", fontweight="bold")
ax.text(2.9, 15.3, "Gaussian 0.45", color=MAG, fontsize=8.8, ha="left", va="bottom", fontweight="bold")
ax.text(16.3, 18.55, "homogeneous 0.45", color=MAG, fontsize=8.8, ha="left", va="top", fontweight="bold")
ax.text(17.2, 16.4, "experiment", color=INK, fontsize=8.8, ha="left", va="top", fontweight="bold")
ax.set_xlim(0, 24.2); ax.set_ylim(0, 25.5)
ax.set_xticks(range(0, 21, 2)); ax.set_xlabel("pulse number", color=INK)
ax.set_ylabel("crack extent from the borehole wall (mm)", color=INK)
ax.set_title(f"AT1, spectral split, absorbing side faces. Solid lines/circles: Gaussian Gc scatter (cv 0.12); dashed/triangles: homogeneous Gc; colour = eta.\n"
             f"Filled marker: end of a completed 60 us window; hollow: pulse in progress. Read {datetime.datetime.now():%Y-%m-%d, %H:%M} PDT.", color=INK2, fontsize=7.8, loc="left")
fig.suptitle("Radial crack growth per pulse: every 2D section case against the acrylic-cube measurements", x=0.02, ha="left", color=INK, fontsize=11)
ax.grid(True, color=GRID, lw=0.8); ax.set_axisbelow(True)
for s in ("top", "right"): ax.spines[s].set_visible(False)
for s in ("left", "bottom"): ax.spines[s].set_color(GRID)
ax.tick_params(colors=INK2, labelsize=9)
ax.legend(handles, labels, loc="lower right", frameon=False, fontsize=7.8)
fig.tight_layout(); fig.savefig(OUT, facecolor=PAPER); print("wrote", OUT)
