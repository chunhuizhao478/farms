# Paper figure: the borehole pressure history of the model. (a) one pulse, the double-exponential shape with the peak P_p at the
# rise time t_r and 10 % of the peak at t_d; (b) the train of 60 us windows with the damping schedule (zeta_pulse during the
# first 25 us of each window, the quench value afterwards).  python3 fig_pulse_history.py OUT_STEM  -> OUT_STEM.pdf and .png
import sys, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
STEM = sys.argv[1]
plt.rcParams.update({"font.size": 9, "mathtext.fontset": "dejavusans", "axes.linewidth": 0.8, "pdf.fonttype": 42, "ps.fonttype": 42})
a, b, tr, td, Tw, tfree = 4.1032e5, 1.1529e7, 0.3e-6, 6e-6, 60e-6, 25e-6   # 1/s, 1/s, s, s, s, s
Pp = 214.4                                                                 # MPa, eta = 0.45 (p_wall_max of the run)
def p(t): return Pp * (np.exp(-a * t) - np.exp(-b * t)) / (np.exp(-a * tr) - np.exp(-b * tr))
BLUE, INK, GREY = "#2a78d6", "#000000", "#555555"
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(5.2, 5.0), dpi=300, gridspec_kw={"height_ratios": [1, 1], "hspace": 0.55})
# (a) one pulse
t = np.linspace(0, 12e-6, 4000)
ax1.plot(t * 1e6, p(t), color=BLUE, lw=1.6)
for x, lab in ((tr, "$t_r$"), (td, "$t_d$")):
    ax1.axvline(x * 1e6, color=GREY, lw=0.8, ls=(0, (3, 3)))
    ax1.text(x * 1e6 + 0.12, Pp * 1.02, lab, color=INK, fontsize=9, ha="left", va="bottom")   # values of t_r, t_d, P_p, a, b: caption or table
ax1.axhline(0.1 * Pp, color=GREY, lw=0.8, ls=(0, (1, 2)))
ax1.text(11.9, 0.1 * Pp + 4, "$0.1\\,P_p$", color=INK, fontsize=8.5, ha="right", va="bottom")
ax1.text(11.9, Pp * 0.98, "$p(t) = P_p\\,\\dfrac{e^{-at} - e^{-bt}}{e^{-a t_r} - e^{-b t_r}}$", color=INK, fontsize=9, ha="right", va="top")
ax1.set_xlim(0, 12); ax1.set_ylim(0, 240); ax1.set_yticks(range(0, 241, 60))
ax1.set_xlabel("Time (µs)"); ax1.set_ylabel("Borehole pressure (MPa)")
ax1.set_title("(a) one pulse", loc="left", fontsize=9.5)
# (b) the train
t = np.linspace(0, 3 * Tw, 30000); pt = np.zeros_like(t)
for k in range(3):
    m = t >= k * Tw; pt[m] = p(t[m] - k * Tw)
for k in range(3):
    ax2.axvspan(k * Tw * 1e6, (k * Tw + tfree) * 1e6, color="#dbe7f7", lw=0, zorder=0)
    ax2.axvspan((k * Tw + tfree) * 1e6, (k + 1) * Tw * 1e6, color="#ebebe8", lw=0, zorder=0)
    ax2.text((k * Tw + 3e-6) * 1e6, Pp * 0.92, f"pulse {k + 1}", color=INK, fontsize=8.5, ha="left", va="top")
ax2.plot(t * 1e6, pt, color=BLUE, lw=1.3, zorder=3)
ax2.set_xlim(0, 180); ax2.set_ylim(0, 240); ax2.set_yticks(range(0, 241, 60)); ax2.set_xticks(range(0, 181, 30))
ax2.set_xlabel("Time (µs)"); ax2.set_ylabel("Borehole pressure (MPa)")
ax2.set_title("(b) pulse train, one 60 µs window per pulse", loc="left", fontsize=9.5)
ax2.text(15.5, 100, "$\\zeta = 10^{-8}$ s", color=INK, fontsize=8.5, ha="center", va="center")
ax2.text(42.5, 100, "$\\zeta = 10^{-5}$ s\n(quench)", color=INK, fontsize=8.5, ha="center", va="center", linespacing=1.4)
for ax in (ax1, ax2):
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    ax.tick_params(labelsize=8.5, length=3)
fig.savefig(STEM + ".pdf", bbox_inches="tight", pad_inches=0.03); fig.savefig(STEM + ".png", dpi=300, bbox_inches="tight", pad_inches=0.03)
print("wrote", STEM + ".pdf", STEM + ".png")
