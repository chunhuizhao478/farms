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

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
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
A = 6_000.0  # (m) depth where the hydrostatic-to-lithostatic transition *begins*
B = 8_000.0  # (m) depth where the hydrostatic-to-lithostatic transition *ends*
C = 7_500.0  # (m) depth where super-lithostatic transition *begins* (Pf→Sv starts)
D = 12_000.0  # (m) depth where super-lithostatic transition *ends* (Pf≈Sv)
z_max = 20_000.0  # (m) modelling depth
n_pts = 600  # number of depth samples

density_fluid = 1_000.0  # (kg/m³) pore-fluid density
rho = 2_670.0  # (kg/m³) bulk rock density (lithostatic)
g = 9.8  # (m/s²) gravitational acceleration

# Stress coefficients (unchanged)
# b_xx = 0.4
# b_yy = 1.073206
# b_xy = -0.8

# Coefficients for horizontal and shear stresses
b_xx = 0.926793  # 0.4
b_yy = 1.073206
b_xy = -0.8  # -0.8

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
# VERTICAL STRESS (compression negative) - needed for pore pressure calculation
# -------------------------------------------------------------------
sigma_zz = -rho * g * depths
sigma_v_mag = np.abs(sigma_zz)  # magnitude for ratio calculations

# -------------------------------------------------------------------
# PORE PRESSURE (piecewise with super-lithostatic transition)
# -------------------------------------------------------------------
delta_rho = rho - density_fluid  # density contrast
Pf = np.empty_like(depths)

# Region 1: Hydrostatic above A (0 to 6 km)
mask1 = depths <= A
Pf[mask1] = density_fluid * g * depths[mask1]

# Region 2: Hydrostatic-to-lithostatic transition A < z ≤ C (6 km to 7.5 km)
# This region now only goes up to C instead of B
mask2 = (depths > A) & (depths <= C)
z2 = depths[mask2]
Pf_A = density_fluid * g * A
Pf[mask2] = Pf_A + g * (
    density_fluid * (z2 - A) + 0.5 * delta_rho * (z2 - A) ** 2 / (B - A)
)

# Calculate Pf and gradient at depth C (end of Region 2) for continuity
Pf_C = Pf_A + g * (density_fluid * (C - A) + 0.5 * delta_rho * (C - A) ** 2 / (B - A))
# Gradient at C from Region 2: dPf/dz = ρ_fluid·g + δρ·g·(C-A)/(B-A)
grad_C = density_fluid * g + delta_rho * g * (C - A) / (B - A)

# Calculate target Pf and gradient at depth D
lambda_max = 0.95  # target ratio (Pf approaches 95% of Sv)
Pf_D = lambda_max * rho * g * D
grad_D = lambda_max * rho * g  # gradient approaches ρ·g at D

# Region 3: Super-lithostatic transition C < z ≤ D (7.5 km to 12 km)
# Use cubic Hermite interpolation for C1 continuity (smooth gradient)
mask3 = (depths > C) & (depths <= D)
z3 = depths[mask3]
# Normalized coordinate t ∈ [0, 1]
t = (z3 - C) / (D - C)
# Hermite basis functions
h00 = 2 * t**3 - 3 * t**2 + 1  # value at C
h10 = t**3 - 2 * t**2 + t  # gradient at C
h01 = -2 * t**3 + 3 * t**2  # value at D
h11 = t**3 - t**2  # gradient at D
# Cubic Hermite spline
Pf[mask3] = h00 * Pf_C + h10 * (D - C) * grad_C + h01 * Pf_D + h11 * (D - C) * grad_D

# Region 4: Near-lithostatic beyond D (> 12 km)
mask4 = depths > D
z4 = depths[mask4]
Pf[mask4] = lambda_max * rho * g * z4

# -------------------------------------------------------------------
# Print pore pressure ratio diagnostics
# -------------------------------------------------------------------
print("=" * 70)
print("PORE PRESSURE RATIO (λ = Pf / |σv|) AT KEY DEPTHS:")
print("=" * 70)
for z_check in [6000, 7500, 8000, 10000, 12000, 15000, 20000]:
    if z_check <= z_max:
        idx_check = np.argmin(np.abs(depths - z_check))
        lambda_check = Pf[idx_check] / sigma_v_mag[idx_check]
        print(
            f"Depth = {z_check / 1e3:5.1f} km: λ = {lambda_check:.4f}, "
            f"Pf = {Pf[idx_check] / 1e6:6.2f} MPa, |σv| = {sigma_v_mag[idx_check] / 1e6:6.2f} MPa"
        )
print("=" * 70)
print()

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
# SEISMIC PROPERTIES (constant values)
# -------------------------------------------------------------------
vs = 3464 * np.ones_like(depths)  # m/s
vp = 6000 * np.ones_like(depths)  # m/s
rho_depth = 2670 * np.ones_like(depths)  # kg/m^3

# -------------------------------------------------------------------
# PLOT: Seismic Properties vs Depth
# -------------------------------------------------------------------
plt.figure(figsize=(6, 8))
(line_vs,) = plt.plot(vs / 1e3, depths / 1e3, label="Vs (km/s)")
(line_vp,) = plt.plot(vp / 1e3, depths / 1e3, label="Vp (km/s)")
(line_rho,) = plt.plot(
    rho_depth / 1e3, depths / 1e3, label="Density (g/cc)", linestyle="--"
)  # 1 g/cc = 1000 kg/m^3

# Add text labels near the lines with their values
vs_val = float(np.mean(vs) / 1e3)
vp_val = float(np.mean(vp) / 1e3)
rho_val = float(np.mean(rho_depth) / 1e3)
x_off = 0.05  # small horizontal offset so text doesn't sit on the line

# Choose distinct depths (in km) for annotations to avoid overlap
y_vs, y_vp, y_rho = 2.0, 5.0, 8.0

plt.text(
    vs_val + x_off,
    y_vs,
    f"Vs = {vs_val:.2f} km/s",
    color=line_vs.get_color(),
    va="center",
    ha="center",
    rotation=90,
    fontsize=12,
    bbox=dict(facecolor=THEME_LT1, alpha=0.6, edgecolor="none"),
)
plt.text(
    vp_val + x_off,
    y_vp,
    f"Vp = {vp_val:.2f} km/s",
    color=line_vp.get_color(),
    va="center",
    ha="center",
    rotation=90,
    fontsize=12,
    bbox=dict(facecolor=THEME_LT1, alpha=0.6, edgecolor="none"),
)
plt.text(
    rho_val + x_off,
    y_rho,
    f"ρ = {rho_val:.2f} g/cc",
    color=line_rho.get_color(),
    va="center",
    ha="center",
    rotation=90,
    fontsize=12,
    bbox=dict(facecolor=THEME_LT1, alpha=0.6, edgecolor="none"),
)
plt.gca().invert_yaxis()
plt.ylabel("Depth (km)", fontsize=20)
plt.xlabel("Value (km/s or g/cc)", fontsize=20)
plt.title("Seismic Properties vs Depth", fontsize=20)
plt.legend(loc="best", fontsize=18)
ax = plt.gca()
tick_prop = mpl.font_manager.FontProperties(
    family="DejaVu Sans", size=16
)  # e.g., 'Arial', 'Helvetica'
for lab in ax.get_xticklabels() + ax.get_yticklabels():
    lab.set_fontproperties(tick_prop)
plt.tight_layout()
plt.savefig("seismic_properties_vs_depth.png", dpi=300)
plt.show()

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
plt.title("Stress Components vs Depth", fontsize=20)
plt.legend(loc="best", fontsize=10)
# plt.grid(True, which='both', linestyle=':')
# Set tick label font family and size
ax = plt.gca()
tick_prop = mpl.font_manager.FontProperties(
    family="DejaVu Sans", size=16
)  # e.g., 'Arial', 'Helvetica'
for lab in ax.get_xticklabels() + ax.get_yticklabels():
    lab.set_fontproperties(tick_prop)
plt.tight_layout()
plt.savefig("stress_components_vs_depth.png", dpi=300)
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

# -------------------------------------------------------------------
# PLOT: Strain Invariants
# -------------------------------------------------------------------
plt.figure(figsize=(6, 8))
plt.plot(I1, depths, label=r"$I_1$")
plt.plot(I2, depths, label=r"$I_2$")
plt.plot(xi, depths, label=r"$\xi$")
plt.gca().invert_yaxis()
plt.ylabel("Depth (m)")
plt.xlabel("Invariant values")
plt.title("Strain Invariants vs Depth")
plt.legend(loc="best")
plt.grid(True)
plt.tight_layout()
plt.savefig("strain_invariants_vs_depth.png", dpi=300)
plt.show()

# -------------------------------------------------------------------
# PRINCIPAL STRESSES at target depth
# -------------------------------------------------------------------
depth_target = 7500
idx = int(np.argmin(np.abs(depths - depth_target)))
depth_exact = depths[idx]

sxx = sigma_xx[idx]
syy = sigma_yy[idx]
szz = sigma_zz[idx]
sxy = sigma_xy[idx]
Pf_target = Pf[idx]
xi_target = xi[idx]

s_tensor = np.array([[sxx, sxy], [sxy, syy]])
eigvals, eigvecs = np.linalg.eigh(s_tensor)
s2, s1 = eigvals
θ_deg = np.degrees(np.arctan2(eigvecs[1, 0], eigvecs[0, 0]))

print(f"Depth (array snap-to):  {depth_exact / 1e3:.3f} km")
print(f"Pore-fluid pressure:   {Pf_target / 1e6:8.2f}  MPa")
print(f"szz (vertical):        {szz / 1e6:8.2f}  MPa")
print(
    f"sxx, syy, sxy:         {sxx / 1e6:8.2f}, {syy / 1e6:8.2f}, {sxy / 1e6:8.2f} MPa"
)
print(f"s1 (major):            {s1 / 1e6:8.2f}  MPa")
print(f"s2 (minor):            {s2 / 1e6:8.2f}  MPa")
print(f"Orientation of s1:     {θ_deg:6.2f}°  (CCW from +x)")
print(f"Invariant xi:          {xi_target:.3f}")

# -------------------------------------------------------------------
# Compute Principal Stresses along entire depth profile
# For strike-slip system: Sv (vertical), S_Hmax, S_hmin
# -------------------------------------------------------------------
# Initialize arrays for principal stresses
S_Hmax = np.zeros_like(depths)
S_hmin = np.zeros_like(depths)
Sv = sigma_zz  # Vertical/overburden stress

# Compute horizontal principal stresses at each depth
for i in range(len(depths)):
    # 2D stress tensor in horizontal plane (xx-yy plane)
    s_tensor_2d = np.array([[sigma_xx[i], sigma_xy[i]], [sigma_xy[i], sigma_yy[i]]])

    # Compute eigenvalues (principal stresses in horizontal plane)
    eigvals = np.linalg.eigvalsh(s_tensor_2d)  # ascending order

    # For strike-slip: S_Hmax > Sv > S_hmin (in compression convention)
    # Compression is negative, so eigvals[0] is most negative (most compressive)
    S_Hmax[i] = eigvals[0]  # most compressive (most negative, larger magnitude)
    S_hmin[i] = eigvals[1]  # less compressive (less negative, smaller magnitude)

# ------------------------
# Plot Principal Stresses - Combined Figure
# ------------------------
plt.figure(figsize=(6, 8))
plt.plot(
    np.abs(S_Hmax) / 1e6,
    depths / 1e3,
    linewidth=2,
    label="S$_{Hmax}$",
    color=THEME_ACCENTS[0],
)
plt.plot(
    np.abs(Sv) / 1e6, depths / 1e3, linewidth=2, label="S$_v$", color=THEME_ACCENTS[2]
)
plt.plot(
    np.abs(S_hmin) / 1e6,
    depths / 1e3,
    linewidth=2,
    label="S$_{hmin}$",
    color=THEME_ACCENTS[1],
)
plt.plot(
    Pf / 1e6,
    depths / 1e3,
    linewidth=2,
    label="P$_p$",
    linestyle="--",
    color=THEME_ACCENTS[3],
)
plt.gca().invert_yaxis()
plt.ylabel("Depth (km)", fontsize=20)
plt.xlabel("Stress (MPa)", fontsize=20)
plt.title("Principal Stresses vs Depth", fontsize=20)
plt.legend(loc="best", fontsize=18)
plt.grid(True, which="both", linestyle=":")
ax = plt.gca()
tick_prop = mpl.font_manager.FontProperties(family="DejaVu Sans", size=16)
for lab in ax.get_xticklabels() + ax.get_yticklabels():
    lab.set_fontproperties(tick_prop)
plt.tight_layout()
plt.savefig("principal_stresses_combined_vs_depth.png", dpi=300)
plt.show()
