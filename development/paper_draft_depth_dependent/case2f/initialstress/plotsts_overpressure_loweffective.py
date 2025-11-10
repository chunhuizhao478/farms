"""Enhanced pore-pressure and stress profile with a depth-dependent transition.

Between depths A and B (m), the vertical pore-pressure gradient switches smoothly
from hydrostatic (density_fluid · g) to lithostatic/over-pressure (rho · g).
Below B the over-pressured gradient is retained.  All other stress components
follow the same empirical relationships defined in the original *plotsts.py*:

    σ_xx = b_xx (σ_zz + P_f) − P_f
    σ_yy = b_yy (σ_zz + P_f) − P_f
    σ_xy = b_xy (σ_zz + P_f)

with coefficients valid down to 15.6 km.  Cohesion, static, and residual
strengths are computed exactly as before.

Edit the constants under “USER-PARAMETERS” to suit your model.
"""

from __future__ import annotations

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from cycler import cycler

# -------------------------------------------------------------------
# THEME COLORS (from your PPT deck)
# -------------------------------------------------------------------
THEME_ACCENTS = ["#000000", "#577399", "#C22B48", "#FFAB40", "#579499", "#5C5799"]
THEME_LT1 = "#FFFFFF"  # light background
THEME_DK1 = "#000000"  # dark text

mpl.rcParams.update(
    {
        "axes.prop_cycle": cycler("color", THEME_ACCENTS),
        "figure.facecolor": THEME_LT1,
        "axes.facecolor": THEME_LT1,
        "savefig.facecolor": THEME_LT1,
        "text.color": THEME_DK1,
        "axes.labelcolor": THEME_DK1,
        "xtick.color": THEME_DK1,
        "ytick.color": THEME_DK1,
        "axes.edgecolor": THEME_DK1,
        "grid.color": THEME_DK1,
        "legend.edgecolor": THEME_DK1,
    }
)

# -------------------------------------------------------------------
# USER-PARAMETERS
# -------------------------------------------------------------------
A = 6_000.0  # (m) depth where the transition *begins*
B = 8_000.0  # (m) depth where the transition *ends*
z_max = 20_000.0  # (m) modelling depth
n_pts = 600  # number of depth samples

density_fluid = 1_000.0  # (kg/m³) pore-fluid density
rho = 2_670.0  # (kg/m³) bulk rock density (lithostatic)
g = 9.8  # (m/s²) gravitational acceleration

# Pore pressure ratio for depth > B
# Controls how much overburden is supported by pore fluid:
#   0.95 = 95% supported by fluid, 5% effective stress (low effective stress)
#   0.98 = 98% supported by fluid, 2% effective stress (very low)
#   0.99 = 99% supported by fluid, 1% effective stress (extremely low)
lambda_pp = 0.98  # Pore pressure ratio below 8 km (HIGH for low effective stress)

# Stress coefficients (unchanged)
b_xx = 0.926793  # 0.4
b_yy = 1.073206
b_xy = -0.8

# Friction & cohesion parameters (unchanged)
mu_s = 0.8  # static friction coefficient
mu_d = 0.6  # dynamic / residual friction coefficient
c0 = 0.4e6  # Pa, cohesion at 0–4 km
dc_dz = 0.00072e6  # Pa/m, cohesion gradient up to 4 km

# -------------------------------------------------------------------
# DEPTH GRID
# -------------------------------------------------------------------
depths = np.linspace(0.0, z_max, n_pts)  # positive downward (m)

# -------------------------------------------------------------------
# PORE PRESSURE (piecewise-analytic)
# -------------------------------------------------------------------
delta_rho = rho - density_fluid  # density contrast
Pf = np.empty_like(depths)

# Region 1: Hydrostatic above A
mask1 = depths <= A
Pf[mask1] = density_fluid * g * depths[mask1]

# Region 2: Quadratic transition A < z ≤ B
# Transition from hydrostatic (Pf_A) to lambda_pp * overburden (Pf_B_target)
mask2 = (depths > A) & (depths <= B)
z2 = depths[mask2]
Pf_A = density_fluid * g * A
Pf_B_target = lambda_pp * rho * g * B  # Target pore pressure at B
# Quadratic interpolation from Pf_A to Pf_B_target (gradual then steeper)
s = (z2 - A) / (B - A)  # Normalized depth parameter [0, 1]
Pf[mask2] = Pf_A + (Pf_B_target - Pf_A) * s**2

# Region 3: Over-pressured below B (lambda_pp fraction of overburden)
mask3 = depths > B
z3 = depths[mask3]
Pf[mask3] = lambda_pp * rho * g * z3

# -------------------------------------------------------------------
# VERTICAL STRESS (compression negative)
# -------------------------------------------------------------------
sigma_zz = -rho * g * depths

# -------------------------------------------------------------------
# HORIZONTAL & SHEAR STRESSES
# -------------------------------------------------------------------
Omega = np.ones_like(depths)
Omega[(depths > 15000) & (depths <= 20000)] = (
    20000 - depths[(depths > 15000) & (depths <= 20000)]
) / 5000
Omega[depths > 20000] = 0.0

sigma_xx = Omega * (b_xx * (sigma_zz + Pf) - Pf) + (1 - Omega) * sigma_zz
sigma_yy = Omega * (b_yy * (sigma_zz + Pf) - Pf) + (1 - Omega) * sigma_zz
sigma_xy = Omega * (b_xy * (sigma_zz + Pf))

# -------------------------------------------------------------------
# COHESION
# -------------------------------------------------------------------
mask_cohesion = depths <= 5_000.0
c = np.where(mask_cohesion, c0 + dc_dz * (5_000.0 - depths), c0)

# -------------------------------------------------------------------
# SHEAR STRENGTHS
# -------------------------------------------------------------------
static_shear_strength = c + np.abs(mu_s * (sigma_yy + Pf))
residual_shear_strength = c + np.abs(mu_d * (sigma_yy + Pf))

# -------------------------------------------------------------------
# PLOT: Stress Components
# -------------------------------------------------------------------
plt.figure(figsize=(6, 8))
plt.plot(Pf / 1e6, depths / 1e3, label=r"$P_f$")
plt.plot(np.abs(sigma_zz) / 1e6, depths / 1e3, label=r"$\sigma_{zz}$")
plt.plot(np.abs(sigma_xx) / 1e6, depths / 1e3, label=r"$\sigma_{xx}$")
plt.plot(np.abs(sigma_yy) / 1e6, depths / 1e3, label=r"$\sigma_{yy}$")
plt.plot(sigma_xy / 1e6, depths / 1e3, label=r"$\sigma_{xy}$")
plt.plot(
    static_shear_strength / 1e6,
    depths / 1e3,
    label="Static Shear Strength",
    linestyle="--",
    alpha=0.8,
)
plt.plot(
    residual_shear_strength / 1e6,
    depths / 1e3,
    label="Residual Shear Strength",
    linestyle="--",
    alpha=0.8,
)

plt.gca().invert_yaxis()
plt.ylabel("Depth (km)", fontsize=20)
plt.xlabel("Stress (MPa)", fontsize=20)
plt.title("Total Stress Components vs Depth", fontsize=20)
plt.legend(loc="best", fontsize=18)
# plt.grid(True, which='both', linestyle=':')
# Set tick label font family and size
ax = plt.gca()
tick_prop = mpl.font_manager.FontProperties(
    family="DejaVu Sans", size=16
)  # e.g., 'Arial', 'Helvetica'
for lab in ax.get_xticklabels() + ax.get_yticklabels():
    lab.set_fontproperties(tick_prop)
plt.tight_layout()
plt.savefig("total_stress_components_vs_depth.png", dpi=300)
plt.show()

# -------------------------------------------------------------------
# Effective stress plot
sigma_xx = sigma_xx + Pf
sigma_yy = sigma_yy + Pf
sigma_zz = sigma_zz + Pf

static_shear_strength = c + np.abs(mu_s * (sigma_yy))
residual_shear_strength = c + np.abs(mu_d * (sigma_yy))

# -------------------------------------------------------------------
# PLOT: Stress Components
# -------------------------------------------------------------------
plt.figure(figsize=(6, 8))
# plt.plot(Pf / 1e6, depths / 1e3, label=r"$P_f$")
plt.plot(np.abs(sigma_zz) / 1e6, depths / 1e3, label=r"$\sigma_{zz}$")
plt.plot(np.abs(sigma_xx) / 1e6, depths / 1e3, label=r"$\sigma_{xx}$")
plt.plot(np.abs(sigma_yy) / 1e6, depths / 1e3, label=r"$\sigma_{yy}$")
plt.plot(sigma_xy / 1e6, depths / 1e3, label=r"$\sigma_{xy}$")
plt.plot(
    static_shear_strength / 1e6,
    depths / 1e3,
    label="Static Shear Strength",
    linestyle="--",
    alpha=0.8,
)
plt.plot(
    residual_shear_strength / 1e6,
    depths / 1e3,
    label="Residual Shear Strength",
    linestyle="--",
    alpha=0.8,
)

plt.gca().invert_yaxis()
plt.ylabel("Depth (km)", fontsize=20)
plt.xlabel("Stress (MPa)", fontsize=20)
plt.title("Stress Components vs Depth", fontsize=20)
plt.legend(loc="best", fontsize=18)
# plt.grid(True, which='both', linestyle=':')
# Set tick label font family and size
ax = plt.gca()
tick_prop = mpl.font_manager.FontProperties(
    family="DejaVu Sans", size=16
)  # e.g., 'Arial', 'Helvetica'
for lab in ax.get_xticklabels() + ax.get_yticklabels():
    lab.set_fontproperties(tick_prop)
plt.tight_layout()
plt.savefig("effective_stress_components_vs_depth.png", dpi=300)
plt.show()

# -------------------------------------------------------------------
# STRAIN COMPUTATION (constant μ, λ here)
# -------------------------------------------------------------------
mu = 32.04e9
lmbda = 32.04e9

stress = np.zeros((len(depths), 3, 3))
stress[:, 0, 0] = sigma_xx
stress[:, 1, 1] = sigma_yy
stress[:, 2, 2] = sigma_zz
stress[:, 0, 1] = sigma_xy
stress[:, 1, 0] = sigma_xy

trace_coeff = -lmbda / (2 * mu * (3 * lmbda + 2 * mu))
strain = stress / (2.0 * mu)
enum = np.trace(stress, axis1=1, axis2=2)
for k in range(3):
    strain[:, k, k] += trace_coeff * enum

I1 = np.trace(strain, axis1=1, axis2=2)
I2 = np.einsum("nij,nij->n", strain, strain)
xi = I1 / np.sqrt(I2)

# ------------------------
# Plot Strain Invariants
# ------------------------
plt.figure(figsize=(6, 8))
# plt.plot(I1, depths, label=r'$I_1$')
# plt.plot(I2, depths, label=r'$I_2$')
plt.plot(xi, depths / 1e3, label=r"$\xi$")
plt.gca().invert_yaxis()
plt.ylabel("Depth (km)", fontsize=20)
plt.xlabel("Invariant values", fontsize=20)
plt.title("Strain Invariants vs Depth", fontsize=20)
plt.grid(True, which="both", linestyle=":")
ax2 = plt.gca()
# Fixed ticks for xi as requested (-1.5 to 1.5)
ax2.set_xticks([-1.73, -1.075, -0.5, 0.0])
# Optionally enforce symmetric x-limits if xi range narrower
current_xlim = ax2.get_xlim()
if current_xlim[0] > -1.5 or current_xlim[1] < 1.5:
    ax2.set_xlim(-1.8, 0)
for lab in ax2.get_xticklabels() + ax2.get_yticklabels():
    lab.set_fontproperties(tick_prop)
plt.legend(loc="best", fontsize=14)
plt.tight_layout()
plt.savefig("strain_invariants_vs_depth.png", dpi=300)
plt.show()
