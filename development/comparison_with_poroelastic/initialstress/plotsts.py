import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from cycler import cycler

# ----------------------------------------------
# THEME (mirrors case1f/plotsts.py for consistency)
# ----------------------------------------------
THEME_ACCENTS = ["#000000", "#577399", "#C22B48", "#FFAB40", "#579499", "#5C5799"]
THEME_LT1 = "#FFFFFF"
THEME_DK1 = "#000000"

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

# Physical constants
density_fluid = 1000  # kg/m^3
rho = 2670  # kg/m^3
g = 9.8  # m/s^2

xi_o = -1.0
d_0 = 0
gamma_r = 35e9
# Elastic constants
mu = 32.04e9 + xi_o * d_0 * gamma_r  # shear modulus (Pa)
lmbda = 32.04e9  # Lame's first parameter (Pa)

# Depth array (m)
depths = np.linspace(0, 20000, 600)  # from surface to 20 km

# Pore pressure and vertical stress
Pf = density_fluid * g * depths
sigma_zz = -rho * g * depths

# Coefficients for horizontal and shear stresses
b_xx = 0.4  # 0.4
b_yy = 1.073206
b_xy = -0.8  # -0.8

# Coefficients for horizontal and shear stresses
# b_xx = 0.4
# b_yy = 1.073206
# b_xy = -0.8

# Piecewise definitions
# Define tapering coefficient Omega(depth)
Omega = np.ones_like(depths)
Omega[(depths > 15000) & (depths <= 20000)] = (
    20000 - depths[(depths > 15000) & (depths <= 20000)]
) / 5000
Omega[depths > 20000] = 0.0

# Piecewise definitions with tapering applied to deviatoric stress
sigma_xx = Omega * (b_xx * (sigma_zz + Pf) - Pf) + (1 - Omega) * sigma_zz
sigma_yy = Omega * (b_yy * (sigma_zz + Pf) - Pf) + (1 - Omega) * sigma_zz
sigma_xy = Omega * (b_xy * (sigma_zz + Pf))

# cohesion
mask = depths <= 4000
c = np.where(mask, 0.4e6 + (0.00072e6) * (5000 - depths), 0.4e6)  # cohesion in Pa

# shear strength
mu_s = 0.8
static_shear_strength = c + abs(mu_s * (sigma_yy + Pf))
mu_d = 0.6
residual_shear_strength = c + abs(mu_d * (sigma_yy + Pf))

# constant material properties
vs = 3464 * np.ones_like(depths)  # m/s
vp = 6000 * np.ones_like(depths)  # m/s
rho_depth = 2670 * np.ones_like(depths)  # kg/m^3
# ------------------------
# Plot depth variation of Vs, Vp, and density
# ------------------------
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

# ------------------------
# Plot Stress Components
# ------------------------
# plt.figure(figsize=(6, 8))
# plt.plot(Pf/1e6, depths/1e3, label=r'$P_f$')
# plt.plot(abs(sigma_zz)/1e6, depths/1e3, label=r'$\sigma_{zz}$') #when plotting, we use positive values for compression
# plt.plot(abs(sigma_xx)/1e6, depths/1e3, label=r'$\sigma_{xx}$')
# plt.plot(abs(sigma_yy)/1e6, depths/1e3, label=r'$\sigma_{yy}$')
# plt.plot(sigma_xy/1e6, depths/1e3, label=r'$\sigma_{xy}$')
# plt.plot(static_shear_strength/1e6, depths/1e3, label='Static Shear Strength', linestyle='--', color='orange')
# plt.plot(residual_shear_strength/1e6, depths/1e3, label='Residual Shear Strength', linestyle='--', color='red')
# plt.gca().invert_yaxis()
# plt.ylabel('Depth (km)')
# plt.xlabel('Stress (MPa)')
# plt.title('Stress Components vs Depth')
# plt.legend(loc='best')
# plt.grid(True)
# plt.tight_layout()
# plt.savefig('stress_components_vs_depth.png', dpi=300)
# plt.show()

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
plt.legend(loc="best", fontsize=18)
plt.grid(True, which="both", linestyle=":")
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

# ------------------------
# Compute Strains
# ------------------------
# Assemble stress tensor components
stress = np.zeros((len(depths), 3, 3))
stress[:, 0, 0] = sigma_xx
stress[:, 1, 1] = sigma_yy
stress[:, 2, 2] = sigma_zz
stress[:, 0, 1] = sigma_xy
stress[:, 1, 0] = sigma_xy

# Inverse Hooke's law: epsilon_ij = 1/(2μ) s_ij - λ/(2μ(3λ+2μ)) s_kk δ_ij
trace_coeff = -lmbda / (2 * mu * (3 * lmbda + 2 * mu))
strain = np.zeros_like(stress)
enum = np.trace(stress, axis1=1, axis2=2)  # vector of s_kk
for i in range(3):
    for j in range(3):
        strain[:, i, j] = stress[:, i, j] / (2 * mu)
for k in range(3):
    strain[:, k, k] += trace_coeff * enum

# Compute invariants: I1, I2, xi
I1 = np.trace(strain, axis1=1, axis2=2)
I2 = np.einsum("nij,nij->n", strain, strain)
xi = I1 / np.sqrt(I2)

# ------------------------
# Plot Strain Invariants
# ------------------------
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

# -----------------------------------------------------------
# Principal stresses and maximum-principal orientation (2-D)
# -----------------------------------------------------------
# --- pick the target depth (m) ---------------------------------------------
depth_target = 7500  # 7.5 km

# ---------------------------------------------------------------------------
# Find the row in `depths` that is closest to the target
idx = int(np.argmin(np.abs(depths - depth_target)))
depth_exact = depths[idx]  # the exact depth value in your array

# ----- stresses already computed in your script ----------------------------
sxx = sigma_xx[idx]
syy = sigma_yy[idx]
szz = sigma_zz[idx]
sxy = sigma_xy[idx]
Pf_single = Pf[idx]
xi_single = xi[idx]  # invariant xi at the target depth

# ----- principal stresses & orientation ------------------------------------
s_tensor = np.array([[sxx, sxy], [sxy, syy]])

eigvals, eigvecs = np.linalg.eigh(s_tensor)  # ascending order
s2, s1 = eigvals  # s1 = major, s2 = minor
θ_deg = np.degrees(
    np.arctan2(
        eigvecs[1, 0],  # angle of s1
        eigvecs[0, 0],
    )
)

# ----- report ---------------------------------------------------------------
print(f"Depth (array snap-to):  {depth_exact / 1e3:.3f} km")
print(f"Pore-fluid pressure:   {Pf_single / 1e6:8.2f}  MPa")
print(f"szz (vertical):        {szz / 1e6:8.2f}  MPa")
print(
    f"sxx, syy, sxy:         {sxx / 1e6:8.2f}, {syy / 1e6:8.2f}, {sxy / 1e6:8.2f} MPa"
)
print(f"s1 (major):            {s1 / 1e6:8.2f}  MPa")
print(f"s2 (minor):            {s2 / 1e6:8.2f}  MPa")
print(f"Orientation of s1:     {θ_deg:6.2f}°  (CCW from +x)")
print(f"Invariant xi:          {xi_single:.3f}")
# ---------------------------------------------------------------------------

# -----------------------------------------------------------
# Compute Principal Stresses along entire depth profile
# For strike-slip system: Sv (vertical), S_Hmax, S_hmin
# -----------------------------------------------------------
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
