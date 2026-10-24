# Paper figure: crack extent at the end of each pulse, homogeneous-Gc cases (eta 0.4 / 0.45 / 0.5) vs the measurements.
#   python3 plot_extent_vs_pulse_homog.py OUT_STEM      -> OUT_STEM.pdf (vector) and OUT_STEM.png (300 dpi)
#   reads vista_csv/att_homog_eta0p*_absorb.csv; no title (the caption carries it), legend only.
#   python3 plot_extent_vs_pulse_homog.py OUT_STEM --with-0p425   -> working chart with the eta 0.425 threshold check added (green)
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
S = "/private/tmp/claude-501/-Users-chunhuizhao-projects-farms-cdms/4d45f359-e932-438b-a688-5f487c854dc5/scratchpad/"
sys.path.insert(0, S); from pulse_table import load, windows
STEM = sys.argv[1]; WITH_0P425 = "--with-0p425" in sys.argv[2:]
plt.rcParams.update({"font.size": 9, "mathtext.fontset": "dejavusans", "axes.linewidth": 0.8, "pdf.fonttype": 42, "ps.fonttype": 42,
                     "xtick.direction": "out", "ytick.direction": "out", "xtick.major.width": 0.8, "ytick.major.width": 0.8})
INK, GREY = "#000000", "#555555"
BLUE, MAG, GOLD, GREEN = "#2a78d6", "#c4379a", "#b8860b", "#1a9c6e"   # validated palette; colour = eta as in the working charts
CASES = [("att_homog_eta0p5_absorb.csv", r"$\eta = 0.5$", BLUE, "o", 4.8),
         ("att_homog_eta0p45_absorb.csv", r"$\eta = 0.45$", MAG, "^", 5.2),
         ("att_homog_eta0p4_absorb.csv", r"$\eta = 0.4$", GOLD, "D", 4.2)]
if WITH_0P425:   # threshold check between 0.4 and 0.45 (job 1093037); working chart only, the paper figure keeps the three calibration values
    CASES.insert(2, ("att_homog_eta0p425_absorb.csv", r"$\eta = 0.425$", GREEN, "v", 5.0))
paper = [(5, 7.97), (11, 13.90), (15, 15.66), (20, 20.45)]
fig, ax = plt.subplots(figsize=(5.2, 3.7), dpi=300)
ax.axhline(21.4, color=GREY, lw=0.8, ls=(0, (2, 3)), zorder=1)
ax.text(0.5, 21.4 + 0.35, "side face (21.4 mm)", color=GREY, fontsize=7.5, va="bottom", ha="left")
for f, lab, col, mk, ms in CASES:
    rows, _ = load(S + "vista_csv/" + f)
    done = [(k, e * 1e3) for k, e, c, t in windows(rows) if c]
    ax.plot([p for p, _ in done], [e for _, e in done], color=col, lw=1.5, marker=mk, ms=ms, mec="white", mew=0.6, label=lab, zorder=4)
    print(f"{lab}: " + ", ".join(f"{p}: {e:.2f}" for p, e in done))
ax.plot([p for p, _ in paper], [e for _, e in paper], color=INK, lw=1.2, marker="s", ms=5.2, mec="white", mew=0.6, label="experiment", zorder=5)
ax.set_xlim(0, 21); ax.set_ylim(0, 25)
ax.set_xticks(range(0, 21, 2)); ax.set_yticks(range(0, 26, 5))
ax.set_xlabel("Pulse number", color=INK); ax.set_ylabel("Crack extent from the borehole wall (mm)", color=INK)
for s in ("top", "right"): ax.spines[s].set_visible(False)
for s in ("left", "bottom"): ax.spines[s].set_color(INK)
ax.tick_params(colors=INK, labelsize=8.5, length=3)
ax.legend(loc="lower right", frameon=False, fontsize=8.5, handlelength=2.6, borderaxespad=0.6)
fig.tight_layout(pad=0.4)
fig.savefig(STEM + ".pdf"); fig.savefig(STEM + ".png", dpi=300); print("wrote", STEM + ".pdf", STEM + ".png")
