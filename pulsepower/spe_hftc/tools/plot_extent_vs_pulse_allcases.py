# Crack extent at the end of each 60 us pulse window for every Vista case, against the paper: python3 ... OUT.png
import csv, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
S = "/private/tmp/claude-501/-Users-chunhuizhao-projects-farms-cdms/4d45f359-e932-438b-a688-5f487c854dc5/scratchpad/vista_csv/"
OUT = sys.argv[1]
TW = 60.0
def series(path, key):
    R = [{k: float(v) for k, v in r.items()} for r in csv.DictReader(open(path))]
    done, partial = [], None
    for k in range(1, 21):
        win = [r for r in R if (k - 1) * TW < r["time"] * 1e6 <= k * TW]
        if not win: break
        if win[-1]["time"] * 1e6 < k * TW - 1.0:
            partial = (k, win[-1][key] * 1e3, win[-1]["time"] * 1e6); break
        done.append((k, win[-1][key] * 1e3))
    return done, partial
INK, INK2, GRID, PAPER = "#0b0b0b", "#52514e", "#e6e5e1", "#fcfcfb"
CASES = [  # file, key, label, color, linestyle, marker
    ("att_eta0p5_absorb.csv", "fracture_extent", "production: Gaussian Gc, eta 0.5, spectral", "#2a78d6", "-", "o"),
    ("att_homog_eta0p5_absorb.csv", "fracture_extent", "homogeneous Gc, eta 0.5", "#c4379a", "-", "^"),
    ("att_homog_eta0p4_absorb.csv", "fracture_extent", "homogeneous Gc, eta 0.4", "#b8860b", "-", "D"),
    ("att_eta0p5_voldev_absorb.csv", "fracture_extent_core", "vol-dev split, Gaussian Gc, eta 0.5 (core measure)", "#52514e", "-.", "x"),
    ("att_eta0p5_free.csv", "fracture_extent_core", "free side faces, Gaussian Gc, eta 0.5 (core measure)", "#15995f", "--", "o"),
]
paper = [(5, 7.97), (11, 13.90), (15, 15.66), (20, 20.45)]
fig, ax = plt.subplots(figsize=(9.2, 5.8), dpi=150); fig.patch.set_facecolor(PAPER); ax.set_facecolor(PAPER)
ax.axhline(21.4, color=INK2, lw=1, ls=(0, (2, 3)), zorder=1)
ax.text(0.4, 21.4 + 0.3, "nearest side face: 21.4 mm from the wall", color=INK2, fontsize=8.5, va="bottom", ha="left")
ax.axhline(10.85, color=INK2, lw=0.8, ls=(0, (1, 3)), zorder=1)
ax.text(23.4, 10.85 + 0.3, "cap of the core measure: 10.85 mm (15 mm radius)", color=INK2, fontsize=8, va="bottom", ha="right")
handles, labels = [], []
for f, key, lab, col, ls, mk in CASES:
    done, part = series(S + f, key)
    if f.startswith("att_eta0p5_free"): done, part = done[:3], None   # the free-face sample fragments by pulse 5
    ms = 7 if mk != "x" else 8
    h, = ax.plot([p for p, _ in done], [e for _, e in done], color=col, lw=2, ls=ls, marker=mk, ms=ms, mec=PAPER if mk != "x" else col, mew=1 if mk != "x" else 2, zorder=4)
    if part:
        k, e, t = part
        ax.plot([done[-1][0] if done else k, k], [done[-1][1] if done else e, e], color=col, lw=1.2, ls=":", zorder=3)
        ax.plot([k], [e], ls="none", marker=mk, ms=ms, mfc=PAPER, mec=col, mew=1.8, zorder=5)
    handles.append(h); labels.append(lab)
    last = part if part else (done[-1][0], done[-1][1], None)
    print(f"{lab}: " + ", ".join(f"{p}: {e:.2f}" for p, e in done) + (f"  | partial pulse {part[0]} at {part[2]:.1f} us: {part[1]:.2f}" if part else ""))
h3, = ax.plot([p for p, _ in paper], [e for _, e in paper], color=INK, lw=1.3, marker="s", ms=7.5, mec=PAPER, mew=1, zorder=6)
for p, e in paper: ax.text(p, e - 0.9, f"{e:.1f}", color=INK2, fontsize=8, ha="center", va="top")
handles.append(h3); labels.append("experiment (Jin et al.): 144 J pulses")
# direct labels
ax.text(17.6, 23.0, "production", color="#2a78d6", fontsize=9, ha="right", va="bottom", fontweight="bold")
ax.annotate("homogeneous eta 0.5", xy=(5.15, 10.0), xytext=(6.6, 7.4), color="#c4379a", fontsize=9, ha="left", va="center", fontweight="bold",
            arrowprops=dict(arrowstyle="-", color="#c4379a", lw=0.8, shrinkA=0, shrinkB=3))
ax.text(6.3, 5.6, "homogeneous eta 0.4", color="#b8860b", fontsize=9, ha="left", va="top", fontweight="bold")
ax.text(19.9, 18.5, "experiment", color=INK, fontsize=9, ha="right", va="top", fontweight="bold")
ax.set_xlim(0, 23.5); ax.set_ylim(0, 25); ax.set_xticks(range(0, 21, 2))
ax.set_xlabel("pulse number", color=INK); ax.set_ylabel("crack extent from the borehole wall (mm)", color=INK)
import datetime
ax.set_title(f"Filled markers: end of a completed 60 us window; hollow: the pulse in progress (or where the run stopped). Read {datetime.datetime.now():%Y-%m-%d, %H:%M} PDT.",
             color=INK2, fontsize=8.3, loc="left")
fig.suptitle("Radial crack growth per pulse: every 2D section case against the acrylic-cube measurements", x=0.02, ha="left", color=INK, fontsize=11)
ax.grid(True, color=GRID, lw=0.8); ax.set_axisbelow(True)
for s in ("top", "right"): ax.spines[s].set_visible(False)
for s in ("left", "bottom"): ax.spines[s].set_color(GRID)
ax.tick_params(colors=INK2, labelsize=9)
ax.legend(handles, labels, loc="lower right", frameon=False, fontsize=8.2)
fig.tight_layout(); fig.savefig(OUT, facecolor=PAPER); print("wrote", OUT)
