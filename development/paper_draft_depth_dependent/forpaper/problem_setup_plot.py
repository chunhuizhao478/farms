"""
Combined problem-setup figure for the paper (2x2 layout).

  (a) Seismic properties          - Vs, Vp, rho
  (b) Background effective stress  - sigma'_xx, sigma'_yy, sigma'_zz, sigma_xy
  (c) On-fault stress              - tau_s, tau_d, tau
  (d) Strain invariant ratio xi
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

# ── Global rc parameters for publication quality ──
mpl.rcParams.update(
    {
        "font.family": "serif",
        "mathtext.fontset": "dejavuserif",
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
        "axes.labelsize": 22,
        "axes.titlesize": 22,
        "xtick.labelsize": 18,
        "ytick.labelsize": 18,
        "legend.fontsize": 14,
        "lines.linewidth": 2.5,
    }
)

# ── Common physical parameters ──
g = 9.8  # m/s^2
density_fluid = 1000.0  # kg/m^3
rho_const = 2670.0  # kg/m^3  (cases 1 & 2)
n_pts = 600
z_max = 20_000.0  # m
depths = np.linspace(0, z_max, n_pts)
depth_km = depths / 1e3

# Stress coefficients (all cases)
b_xx, b_yy, b_xy = 0.926793, 1.073206, -0.8
mu_s_fric, mu_d_fric = 0.8, 0.6

# Tapering Omega (all cases)
Omega = np.ones_like(depths)
taper = (depths > 15_000) & (depths <= 20_000)
Omega[taper] = (20_000 - depths[taper]) / 5_000
Omega[depths > 20_000] = 0.0

# Cohesion: cases 1 & 3 (gradient to 4 km), case 2 (gradient to 5 km)
c_13 = np.where(depths <= 4_000, 0.4e6 + 0.00072e6 * (5_000 - depths), 0.4e6)
c_2 = np.where(depths <= 5_000, 0.4e6 + 0.00072e6 * (5_000 - depths), 0.4e6)

# Elastic constants (cases 1 & 2)
mu_el = 32.04e9
lm_el = 32.04e9
tc_12 = -lm_el / (2 * mu_el * (3 * lm_el + 2 * mu_el))


# ── Helper ──
def compute_stress_strain(sigma_zz, Pf, cohesion, mu_arr, lmbda_arr):
    """Return effective stresses, shear stress, strengths, and strain-invariant xi."""
    sxx = Omega * (b_xx * (sigma_zz + Pf) - Pf) + (1 - Omega) * sigma_zz
    syy = Omega * (b_yy * (sigma_zz + Pf) - Pf) + (1 - Omega) * sigma_zz
    sxy = Omega * (b_xy * (sigma_zz + Pf))

    # Effective stresses
    sxx_e = sxx + Pf
    syy_e = syy + Pf
    szz_e = sigma_zz + Pf

    # Mohr-Coulomb shear strengths (effective normal stress)
    static_str = cohesion + np.abs(mu_s_fric * syy_e)
    resid_str = cohesion + np.abs(mu_d_fric * syy_e)

    # Assemble stress tensor
    stress = np.zeros((n_pts, 3, 3))
    stress[:, 0, 0] = sxx_e
    stress[:, 1, 1] = syy_e
    stress[:, 2, 2] = szz_e
    stress[:, 0, 1] = stress[:, 1, 0] = sxy

    # Inverse Hooke: epsilon_ij = s_ij/(2mu) - lambda/(2mu(3lambda+2mu)) s_kk delta_ij
    if np.ndim(mu_arr) == 0:
        strain = stress / (2.0 * mu_arr)
        tc = -lmbda_arr / (2 * mu_arr * (3 * lmbda_arr + 2 * mu_arr))
    else:
        strain = stress / (2.0 * mu_arr[:, None, None])
        tc = -lmbda_arr / (2 * mu_arr * (3 * lmbda_arr + 2 * mu_arr))

    s_kk = np.trace(stress, axis1=1, axis2=2)
    for k in range(3):
        strain[:, k, k] += tc * s_kk

    I1 = np.trace(strain, axis1=1, axis2=2)
    I2 = np.maximum(np.einsum("nij,nij->n", strain, strain), 1e-18)
    xi = I1 / np.sqrt(I2)

    return sxx_e, syy_e, szz_e, sxy, static_str, resid_str, xi


# ══════════════════════════════════════════════════════
# CASE 1 - constant seismic, hydrostatic Pf
# ══════════════════════════════════════════════════════
vs_c1 = 3464.0 * np.ones_like(depths)
vp_c1 = 6000.0 * np.ones_like(depths)
rho_c1 = 2670.0 * np.ones_like(depths)

Pf_c1 = density_fluid * g * depths
szz_c1 = -rho_const * g * depths

sxx_e_c1, syy_e_c1, szz_e_c1, sxy_c1, ss_c1, rs_c1, xi_c1 = compute_stress_strain(
    szz_c1, Pf_c1, c_13, mu_el, lm_el
)

# ══════════════════════════════════════════════════════
# CASE 2 - constant seismic, overpressured Pf
# ══════════════════════════════════════════════════════
vs_c2, vp_c2, rho_c2 = vs_c1.copy(), vp_c1.copy(), rho_c1.copy()

A_op, B_op = 6_000.0, 8_000.0
delta_rho = rho_const - density_fluid
Pf_c2 = np.empty_like(depths)

m1 = depths <= A_op
Pf_c2[m1] = density_fluid * g * depths[m1]

m2 = (depths > A_op) & (depths <= B_op)
Pf_A = density_fluid * g * A_op
Pf_c2[m2] = Pf_A + g * (
    density_fluid * (depths[m2] - A_op)
    + 0.5 * delta_rho * (depths[m2] - A_op) ** 2 / (B_op - A_op)
)

m3 = depths > B_op
Pf_B = density_fluid * g * B_op + 0.5 * g * delta_rho * (B_op - A_op)
Pf_c2[m3] = Pf_B + rho_const * g * (depths[m3] - B_op)

szz_c2 = -rho_const * g * depths

sxx_e_c2, syy_e_c2, szz_e_c2, sxy_c2, ss_c2, rs_c2, xi_c2 = compute_stress_strain(
    szz_c2, Pf_c2, c_2, mu_el, lm_el
)

# ══════════════════════════════════════════════════════
# CASE 3 - depth-varying seismic (TPV32), hydrostatic Pf
# ══════════════════════════════════════════════════════
depth_knots = np.array(
    [0, 500, 1000, 1600, 2400, 3600, 5000, 9000, 11000, 15000], dtype=float
)
vp_knots = np.array(
    [2200, 3000, 3600, 4400, 4800, 5250, 5500, 5750, 6100, 6300], dtype=float
)
vs_knots = np.array(
    [1050, 1400, 1950, 2500, 2800, 3100, 3250, 3450, 3600, 3700], dtype=float
)
rho_knots = np.array(
    [2200, 2450, 2550, 2600, 2600, 2620, 2650, 2720, 2750, 2900], dtype=float
)

vs_c3 = np.interp(depths, depth_knots, vs_knots, left=vs_knots[0], right=vs_knots[-1])
vp_c3 = np.interp(depths, depth_knots, vp_knots, left=vp_knots[0], right=vp_knots[-1])
rho_c3 = np.interp(
    depths, depth_knots, rho_knots, left=rho_knots[0], right=rho_knots[-1]
)

mu_c3 = rho_c3 * vs_c3**2
lm_c3 = np.maximum(rho_c3 * vp_c3**2 - 2.0 * mu_c3, 0.0)

Pf_c3 = density_fluid * g * depths

# Overburden from depth-varying density
rho_avg = 0.5 * (rho_c3[1:] + rho_c3[:-1])
overburden = np.concatenate(([0.0], np.cumsum(rho_avg * np.diff(depths))))
szz_c3 = -g * overburden

sxx_e_c3, syy_e_c3, szz_e_c3, sxy_c3, ss_c3, rs_c3, xi_c3 = compute_stress_strain(
    szz_c3, Pf_c3, c_13, mu_c3, lm_c3
)

# ══════════════════════════════════════════════════════
# FIGURE — 2×2 panels  (text + arrow labels)
# ══════════════════════════════════════════════════════
#   (a) Seismic properties        — colour = quantity, ls = case group
#   (b) Background effective stress — colour = component, ls = case group
#   (c) On-fault stress            — colour = quantity, ls = case group
#   (d) Strain invariant ratio ξ   — same grouping as seismic
#
# Case groupings
#   (a),(d): solid = Cases 1 & 2,  dashed = Case 3
#   (b),(c): solid = Cases 1 & 3,  dashed = Case 2

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(18, 14))

# Colours — seismic quantities
QVs, QVp, QRho = "#577399", "#C22B48", "#FFAB40"
# Colours — stress components (panel b)
Czz, Cxx, Cyy, Cxy = "#000000", "#577399", "#C22B48", "#FFAB40"
# Colours — on-fault quantities (panel c)
Ctau, Cts, Ctd = "#000000", "#577399", "#C22B48"
# Colours — xi case groups (panel d)
Xi_12, Xi_3 = "#577399", "#C22B48"

akw = dict(arrowstyle="->", lw=1.5)
bbox_kw = dict(facecolor="white", alpha=0.85, edgecolor="none", pad=2)
LZ = 1  # line zorder (lower so annotations sit on top)


def _idx(z_m):
    """Index of nearest depth sample to *z_m* (metres)."""
    return int(np.argmin(np.abs(depths - z_m)))


# ──── (a) Seismic Properties ─────────────────────────
ax1.plot(vs_c1 / 1e3, depth_km, color=QVs, ls="-", lw=2.5, zorder=LZ)
ax1.plot(vs_c3 / 1e3, depth_km, color=QVs, ls="--", lw=2.5, zorder=LZ)
ax1.plot(vp_c1 / 1e3, depth_km, color=QVp, ls="-", lw=2.5, zorder=LZ)
ax1.plot(vp_c3 / 1e3, depth_km, color=QVp, ls="--", lw=2.5, zorder=LZ)
ax1.plot(rho_c1 / 1e3, depth_km, color=QRho, ls="-", lw=2.5, zorder=LZ)
ax1.plot(rho_c3 / 1e3, depth_km, color=QRho, ls="--", lw=2.5, zorder=LZ)

ax1.invert_yaxis()
ax1.set_xlabel(r"Value (km/s or g/cm$^3$)")
ax1.set_ylabel("Depth (km)")
ax1.set_title(r"Seismic Properties ($\rho$, $V_s$, $V_p$)")
ax1.grid(True, which="both", ls=":", alpha=0.5)

ls_h = [Line2D([], [], color="gray", ls="-", lw=2.5),
        Line2D([], [], color="gray", ls="--", lw=2.5)]
ax1.legend(ls_h, ["Cases 1 & 2", "Case 3"],
           loc="lower left", fontsize=14, framealpha=0.9, edgecolor="none")

ax1.annotate(r"$V_s$",
             xy=(vs_c1[_idx(10e3)] / 1e3, 10), xytext=(4.5, 9),
             fontsize=17, fontweight="bold", color=QVs,
             arrowprops=dict(**akw, color=QVs), bbox=bbox_kw)
ax1.annotate(r"$V_p$",
             xy=(vp_c1[_idx(4e3)] / 1e3, 4), xytext=(5.2, 2),
             fontsize=17, fontweight="bold", color=QVp,
             arrowprops=dict(**akw, color=QVp), bbox=bbox_kw)
ax1.annotate(r"$\rho$",
             xy=(rho_c1[_idx(12e3)] / 1e3, 12), xytext=(1.5, 13),
             fontsize=17, fontweight="bold", color=QRho,
             arrowprops=dict(**akw, color=QRho), bbox=bbox_kw)

# ──── (b) Background Effective Stress ─────────────────
ax2.plot(np.abs(szz_e_c1) / 1e6, depth_km, color=Czz, ls="-", lw=2.5, zorder=LZ)
ax2.plot(np.abs(szz_e_c2) / 1e6, depth_km, color=Czz, ls="--", lw=2.5, zorder=LZ)
ax2.plot(np.abs(sxx_e_c1) / 1e6, depth_km, color=Cxx, ls="-", lw=2.5, zorder=LZ)
ax2.plot(np.abs(sxx_e_c2) / 1e6, depth_km, color=Cxx, ls="--", lw=2.5, zorder=LZ)
ax2.plot(np.abs(syy_e_c1) / 1e6, depth_km, color=Cyy, ls="-", lw=2.5, zorder=LZ)
ax2.plot(np.abs(syy_e_c2) / 1e6, depth_km, color=Cyy, ls="--", lw=2.5, zorder=LZ)
ax2.plot(np.abs(sxy_c1) / 1e6, depth_km, color=Cxy, ls="-", lw=2.5, zorder=LZ)
ax2.plot(np.abs(sxy_c2) / 1e6, depth_km, color=Cxy, ls="--", lw=2.5, zorder=LZ)

ax2.invert_yaxis()
ax2.set_xlabel("Stress (MPa)")
ax2.set_ylabel("Depth (km)")
ax2.set_title(r"Background Effective Stress ($\boldsymbol{\sigma}'$)")
ax2.grid(True, which="both", ls=":", alpha=0.5)

ls_h2 = [Line2D([], [], color="gray", ls="-", lw=2.5),
         Line2D([], [], color="gray", ls="--", lw=2.5)]
ax2.legend(ls_h2, ["Cases 1 & 3", "Case 2"],
           loc="upper right", fontsize=14, framealpha=0.9, edgecolor="none")

# Labels in the clear gap between solid lines (right) and dashed Case 2 (≈0)
i8 = _idx(8e3)
i11 = _idx(11e3)
i14 = _idx(14e3)
i10 = _idx(10e3)

ax2.annotate(r"$|\sigma'_{yy}|$",
             xy=(np.abs(syy_e_c1[i8]) / 1e6, 8),
             xytext=(50, 8),
             fontsize=15, fontweight="bold", color=Cyy,
             arrowprops=dict(**akw, color=Cyy), bbox=bbox_kw)
ax2.annotate(r"$|\sigma'_{zz}|$",
             xy=(np.abs(szz_e_c1[i11]) / 1e6, 11),
             xytext=(50, 11),
             fontsize=15, fontweight="bold", color=Czz,
             arrowprops=dict(**akw, color=Czz), bbox=bbox_kw)
ax2.annotate(r"$|\sigma'_{xx}|$",
             xy=(np.abs(sxx_e_c1[i14]) / 1e6, 14),
             xytext=(50, 14),
             fontsize=15, fontweight="bold", color=Cxx,
             arrowprops=dict(**akw, color=Cxx), bbox=bbox_kw)
i16 = _idx(16e3)
ax2.annotate(r"$\sigma_{xy}$",
             xy=(np.abs(sxy_c1[i16]) / 1e6, 16),
             xytext=(250, 18),
             fontsize=15, fontweight="bold", color=Cxy,
             arrowprops=dict(**akw, color=Cxy), bbox=bbox_kw)

# ──── (c) On-Fault Stress ────────────────────────────
ax3.plot(np.abs(sxy_c1) / 1e6, depth_km, color=Ctau, ls="-", lw=2.5, zorder=LZ)
ax3.plot(np.abs(sxy_c2) / 1e6, depth_km, color=Ctau, ls="--", lw=2.5, zorder=LZ)
ax3.plot(ss_c1 / 1e6, depth_km, color=Cts, ls="-", lw=2.5, zorder=LZ)
ax3.plot(ss_c2 / 1e6, depth_km, color=Cts, ls="--", lw=2.5, zorder=LZ)
ax3.plot(rs_c1 / 1e6, depth_km, color=Ctd, ls="-", lw=2.5, zorder=LZ)
ax3.plot(rs_c2 / 1e6, depth_km, color=Ctd, ls="--", lw=2.5, zorder=LZ)

ax3.invert_yaxis()
ax3.set_xlabel("Stress (MPa)")
ax3.set_ylabel("Depth (km)")
ax3.set_title(r"On-Fault Stress ($\tau_s$, $\tau_d$, $\tau$)")
ax3.grid(True, which="both", ls=":", alpha=0.5)

ls_h3 = [Line2D([], [], color="gray", ls="-", lw=2.5),
         Line2D([], [], color="gray", ls="--", lw=2.5)]
ax3.legend(ls_h3, ["Cases 1 & 3", "Case 2"],
           loc="upper right", fontsize=14, framealpha=0.9, edgecolor="none")

# Labels in upper-right clear area (high stress, shallow depth)
i8 = _idx(8e3)
i10 = _idx(10e3)
i12 = _idx(12e3)

ax3.annotate(r"$\tau_s = c + \mu_s |\sigma'_N|$",
             xy=(ss_c1[i8] / 1e6, 8),
             xytext=(175, 2.5),
             fontsize=14, fontweight="bold", color=Cts,
             arrowprops=dict(**akw, color=Cts), bbox=bbox_kw)
ax3.annotate(r"$\tau = |\sigma_{xy}|$",
             xy=(np.abs(sxy_c1[i10]) / 1e6, 10),
             xytext=(175, 5),
             fontsize=15, fontweight="bold", color=Ctau,
             arrowprops=dict(**akw, color=Ctau), bbox=bbox_kw)
ax3.annotate(r"$\tau_d = c + \mu_d |\sigma'_N|$",
             xy=(rs_c1[i12] / 1e6, 12),
             xytext=(175, 7.5),
             fontsize=14, fontweight="bold", color=Ctd,
             arrowprops=dict(**akw, color=Ctd), bbox=bbox_kw)

# Double arrow showing stress drop Δτ between τ (black) and τ_d (red)
i_arr = _idx(14e3)
x_tau = np.abs(sxy_c1[i_arr]) / 1e6   # τ  (black solid)
x_taud = rs_c1[i_arr] / 1e6            # τ_d (red solid)
y_arr = 14.0
ax3.annotate("", xy=(x_tau, y_arr), xytext=(x_taud, y_arr),
             arrowprops=dict(arrowstyle="<->", color=QRho, lw=2),
             zorder=5)
ax3.text(x_taud + (x_tau - x_taud) * 0.15, y_arr - 0.6, r"$\Delta\tau$",
         fontsize=16, ha="center", va="bottom", color=QRho, bbox=bbox_kw)

# ──── (d) Strain Invariant Ratio ξ ───────────────────
ax4.plot(xi_c1, depth_km, color=Xi_12, ls="-", lw=2.5, zorder=LZ)
ax4.plot(xi_c3, depth_km, color=Xi_3, ls="--", lw=2.5, zorder=LZ)

ax4.invert_yaxis()
ax4.set_xlabel(r"$\xi = I_1 / I_2\ (t=0)$")
ax4.set_ylabel("Depth (km)")
ax4.set_title(r"Initial Strain Invariant Ratio ($\xi(t=0)$)")
ax4.set_xlim(-1.8, 0)
ax4.grid(True, which="both", ls=":", alpha=0.5)

ls_h4 = [Line2D([], [], color="gray", ls="-", lw=2.5),
         Line2D([], [], color="gray", ls="--", lw=2.5)]
ax4.legend(ls_h4, ["Cases 1 & 2", "Case 3"],
           loc="lower right", fontsize=14, framealpha=0.9, edgecolor="none")

i8 = _idx(8e3)
i1 = _idx(1e3)
ax4.annotate("Cases 1 & 2",
             xy=(xi_c1[i8], 8), xytext=(-0.4, 10),
             fontsize=15, color=Xi_12,
             arrowprops=dict(**akw, color=Xi_12), bbox=bbox_kw)
ax4.annotate("Case 3",
             xy=(xi_c3[i1], 1), xytext=(-0.4, 3),
             fontsize=15, color=Xi_3,
             arrowprops=dict(**akw, color=Xi_3), bbox=bbox_kw)

# Subplot labels just outside the bottom-left corner of each panel
for ax, label in [(ax1, "(a)"), (ax2, "(b)"), (ax3, "(c)"), (ax4, "(d)")]:
    ax.text(-0.06, -0.01, label, transform=ax.transAxes, fontsize=20,
            va="top", ha="right")

plt.tight_layout()
plt.savefig("problem_setup_combined.png", dpi=300, bbox_inches="tight")
plt.show()
