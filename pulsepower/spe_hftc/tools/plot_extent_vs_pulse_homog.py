# Crack extent at the end of each pulse, homogeneous-Gc cases only (eta 0.4 / 0.45 / 0.5) vs the acrylic-cube measurements:
#   python3 plot_extent_vs_pulse_homog.py OUT.png      (reads vista_csv/att_homog_eta0p*_absorb.csv)
import sys, datetime
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
S = "/private/tmp/claude-501/-Users-chunhuizhao-projects-farms-cdms/4d45f359-e932-438b-a688-5f487c854dc5/scratchpad/"
sys.path.insert(0, S); from pulse_table import load, windows
OUT = sys.argv[1]
INK, INK2, GRID, PAPER = "#0b0b0b", "#52514e", "#e6e5e1", "#fcfcfb"
BLUE, MAG, GOLD = "#2a78d6", "#c4379a", "#b8860b"   # validated palette (dataviz validator, light mode); colour = eta as in the all-cases chart
CASES = [  # file, label, colour, marker
    ("att_homog_eta0p5_absorb.csv", "eta 0.5", BLUE, "o"),
    ("att_homog_eta0p45_absorb.csv", "eta 0.45", MAG, "^"),
    ("att_homog_eta0p4_absorb.csv", "eta 0.4", GOLD, "D"),
]
paper = [(5, 7.97), (11, 13.90), (15, 15.66), (20, 20.45)]
fig, ax = plt.subplots(figsize=(8.6, 5.6), dpi=150)
fig.patch.set_facecolor(PAPER); ax.set_facecolor(PAPER)
ax.axhline(21.4, color=INK2, lw=0.9, ls=(0, (2, 3)), zorder=1)
ax.text(0.4, 21.4 + 0.3, "nearest side face: 21.4 mm from the wall", color=INK2, fontsize=8.5, va="bottom", ha="left")
handles, labels = [], []
for f, lab, col, mk in CASES:
    rows, _ = load(S + "vista_csv/" + f); w = windows(rows)
    done = [(k, e * 1e3) for k, e, c, t in w if c]; part = [(k, e * 1e3, t) for k, e, c, t in w if not c]
    h, = ax.plot([p for p, _ in done], [e for _, e in done], color=col, lw=2, marker=mk, ms=6.5 if mk == "o" else 6, mec=PAPER, mew=1, zorder=4)
    if part:
        k, e, t = part[0]
        ax.plot([done[-1][0], k], [done[-1][1], e], color=col, lw=1.2, ls=(0, (1, 2)), zorder=3)
        ax.plot([k], [e], ls="none", marker=mk, ms=6.5, mfc=PAPER, mec=col, mew=1.6, zorder=6)
        lab += f" (pulse {k} in progress, {t:.0f} us)"
    handles.append(h); labels.append(lab)
    print(f"{lab}: " + ", ".join(f"{p}: {e:.2f}" for p, e in done))
hp, = ax.plot([p for p, _ in paper], [e for _, e in paper], color=INK, lw=1.5, marker="s", ms=7, mec=PAPER, mew=1, zorder=5)
for p, e in paper:   # measured values; the first one sits to the right of its marker, clear of the eta 0.45 curve
    if p == 5: ax.text(p + 0.4, e - 1.0, f"{e:.1f}", color=INK2, fontsize=8, ha="left", va="top")
    else: ax.text(p, e - 0.9, f"{e:.1f}", color=INK2, fontsize=8, ha="center", va="top")
handles.append(hp); labels.append("experiment (Jin et al.): 144 J pulses, 80 s apart")
ax.text(20.35, 23.3, "eta 0.5", color=BLUE, fontsize=9.5, ha="left", va="center", fontweight="bold")
ax.text(14.3, 17.0, "eta 0.45", color=MAG, fontsize=9.5, ha="left", va="center", fontweight="bold")
ax.text(12.5, 9.6, "eta 0.4", color=GOLD, fontsize=9.5, ha="left", va="top", fontweight="bold")
ax.text(17.2, 16.4, "experiment", color=INK, fontsize=9.5, ha="left", va="top", fontweight="bold")
ax.set_xlim(0, 23.5); ax.set_ylim(0, 25.5)
ax.set_xticks(range(0, 21, 2)); ax.set_xlabel("pulse number", color=INK)
ax.set_ylabel("crack extent from the borehole wall (mm)", color=INK)
ax.set_title(f"Homogeneous Gc (2721 J/m2), AT1, spectral split, absorbing side faces.\n"
             f"Extent = farthest d >= 0.5 from the wall at the end of each 60 us window. Read {datetime.datetime.now():%Y-%m-%d, %H:%M} PDT.", color=INK2, fontsize=7.8, loc="left")
fig.suptitle("Radial crack growth per pulse: 2D section model at three eta values against the acrylic-cube measurements", x=0.02, ha="left", color=INK, fontsize=10.5)
ax.grid(True, color=GRID, lw=0.8); ax.set_axisbelow(True)
for s in ("top", "right"): ax.spines[s].set_visible(False)
for s in ("left", "bottom"): ax.spines[s].set_color(GRID)
ax.tick_params(colors=INK2, labelsize=9)
ax.legend(handles, labels, loc="lower right", frameon=False, fontsize=8.5)
fig.tight_layout(); fig.savefig(OUT, facecolor=PAPER); print("wrote", OUT)
