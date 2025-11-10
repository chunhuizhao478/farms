import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from cycler import cycler

# ------------------------
# THEME: PowerPoint palette
# ------------------------
# Accent colors you extracted:
THEME_ACCENTS = ["#000000", "#577399", "#C22B48", "#FFAB40", "#579499", "#5C5799"]

# If you have true lt1/dk1 from your theme XML, put them here.
# Defaults are light background + dark text:
THEME_LT1 = "#FFFFFF"  # light background (figure/axes)
THEME_DK1 = "#000000"  # dark text/ticks/spines

mpl.rcParams.update(
    {
        # line color cycle
        "axes.prop_cycle": cycler("color", THEME_ACCENTS),
        # backgrounds
        "figure.facecolor": THEME_LT1,
        "axes.facecolor": THEME_LT1,
        "savefig.facecolor": THEME_LT1,
        # text / ticks / spines
        "text.color": THEME_DK1,
        "axes.labelcolor": THEME_DK1,
        "xtick.color": THEME_DK1,
        "ytick.color": THEME_DK1,
        "axes.edgecolor": THEME_DK1,
        "grid.color": THEME_DK1,
        "legend.edgecolor": THEME_DK1,
    }
)

# ---------------------------------------------------
# Physics setup
# ---------------------------------------------------
density_fluid = 1000  # kg/m^3
rho = 2670  # kg/m^3 (unused baseline; kept for ref)
g = 9.8  # m/s^2

# Elastic constants (depth-varying via TPV32 profile)
# Table knots (depth in meters)
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

# Depth array (m)
depths = np.linspace(0, 20000, 600)  # from surface to 20 km

# Interpolate seismic properties at all depths (piecewise-linear; clamp at ends)
vp = np.interp(depths, depth_knots, vp_knots, left=vp_knots[0], right=vp_knots[-1])
vs = np.interp(depths, depth_knots, vs_knots, left=vs_knots[0], right=vs_knots[-1])
rho_depth = np.interp(
    depths, depth_knots, rho_knots, left=rho_knots[0], right=rho_knots[-1]
)

# Derived elastic constants (arrays)
mu_depth = rho_depth * vs**2
lmbda_depth = np.maximum(rho_depth * vp**2 - 2.0 * mu_depth, 0.0)

# Pore pressure and vertical stress
Pf = density_fluid * g * depths  # linear with depth for constant fluid density

# Overburden from depth-varying rock density
rho_avg = 0.5 * (rho_depth[1:] + rho_depth[:-1])
dz = np.diff(depths)
overburden = np.concatenate(([0.0], np.cumsum(rho_avg * dz)))  # ∫ rho(z) dz
sigma_zz = -g * overburden  # compressive negative

# Coefficients for horizontal and shear stresses
b_xx = 0.926793
b_yy = 1.073206
b_xy = -0.8

# Tapering coefficient Omega(depth)
Omega = np.ones_like(depths)
mask_taper = (depths > 15000) & (depths <= 20000)
Omega[mask_taper] = (20000 - depths[mask_taper]) / 5000
Omega[depths > 20000] = 0.0

# Deviatoric stress with tapering
sigma_xx = Omega * (b_xx * (sigma_zz + Pf) - Pf) + (1 - Omega) * sigma_zz
sigma_yy = Omega * (b_yy * (sigma_zz + Pf) - Pf) + (1 - Omega) * sigma_zz
sigma_xy = Omega * (b_xy * (sigma_zz + Pf))

# Cohesion (Pa)
mask_shallow = depths <= 4000
c = np.where(mask_shallow, 0.4e6 + (0.00072e6) * (5000 - depths), 0.4e6)

# Shear strengths
mu_s = 0.8
static_shear_strength = c + abs(mu_s * (sigma_yy + Pf))
mu_d = 0.6
residual_shear_strength = c + abs(mu_d * (sigma_yy + Pf))

# ------------------------
# Plot Stress Components
# ------------------------
# plt.figure(figsize=(6, 8))
# plt.plot(Pf/1e6, depths/1e3, label=r'$P_f$')
# plt.plot(abs(sigma_zz)/1e6, depths/1e3, label=r'$\sigma_{zz}$')  # plot compression positive
# plt.plot(abs(sigma_xx)/1e6, depths/1e3, label=r'$\sigma_{xx}$')
# plt.plot(abs(sigma_yy)/1e6, depths/1e3, label=r'$\sigma_{yy}$')
# plt.plot(sigma_xy/1e6, depths/1e3, label=r'$\sigma_{xy}$')
# plt.plot(static_shear_strength/1e6, depths/1e3, label='Static Shear Strength', linestyle='--')
# plt.plot(residual_shear_strength/1e6, depths/1e3, label='Residual Shear Strength', linestyle='--')
# plt.gca().invert_yaxis()
# plt.ylabel('Depth (km)')
# plt.xlabel('Stress (MPa)')
# plt.title('Stress Components vs Depth')
# plt.legend(loc='best', ncol=2)
# plt.grid(True, alpha=0.3)
# plt.tight_layout()
# plt.savefig('stress_components_vs_depth.png', dpi=300)
# plt.show()

# Combined figure with seismic properties, effective stress, and strain invariants side by side
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 8))
tick_prop = mpl.font_manager.FontProperties(family="DejaVu Sans", size=16)

# Left subplot: Seismic Properties (depth-varying for case3f)
ax1.plot(vs / 1e3, depths / 1e3, label="Vs (km/s)")
ax1.plot(vp / 1e3, depths / 1e3, label="Vp (km/s)")
ax1.plot(rho_depth / 1e3, depths / 1e3, label="Density (g/cc)", linestyle="--")
ax1.invert_yaxis()
ax1.set_ylabel("Depth (km)", fontsize=20)
ax1.set_xlabel("Value (km/s or g/cc)", fontsize=20)
ax1.set_title("Seismic Properties vs Depth", fontsize=20)
ax1.legend(loc="best", fontsize=10)
for lab in ax1.get_xticklabels() + ax1.get_yticklabels():
    lab.set_fontproperties(tick_prop)

# ------------------------
# Compute Strains and Effective Stress
# ------------------------

# Convert total stress to effective stress
sigma_xx = sigma_xx + Pf
sigma_yy = sigma_yy + Pf
sigma_zz = sigma_zz + Pf

# Assemble stress tensor components
stress = np.zeros((len(depths), 3, 3))
stress[:, 0, 0] = sigma_xx
stress[:, 1, 1] = sigma_yy
stress[:, 2, 2] = sigma_zz
stress[:, 0, 1] = sigma_xy
stress[:, 1, 0] = sigma_xy

# Inverse Hooke's law with depth-varying μ, λ:
#   epsilon_ij = 1/(2μ) s_ij - λ/(2μ(3λ+2μ)) s_kk δ_ij
trace_coeff = -lmbda_depth / (2 * mu_depth * (3 * lmbda_depth + 2 * mu_depth))
strain = stress / (2.0 * mu_depth[:, None, None])  # first term
enum = np.trace(stress, axis1=1, axis2=2)  # s_kk
for k in range(3):
    strain[:, k, k] += trace_coeff * enum  # add diagonal second term

# Invariants: I1, I2, xi
I1 = np.trace(strain, axis1=1, axis2=2)
I2 = np.einsum("nij,nij->n", strain, strain)
I2 = np.maximum(I2, 1e-18)  # stabilize
xi = I1 / np.sqrt(I2)

# ------------------------
# Third subplot: Strain Invariants (will be added later)
# ------------------------

# # -----------------------------------------------------------
# # Principal stresses and maximum-principal orientation (2-D)
# # -----------------------------------------------------------
# depth_target = 7500  # 7.5 km

# # Find the row in `depths` closest to the target
# idx = int(np.argmin(np.abs(depths - depth_target)))
# depth_exact = depths[idx]

# # Stresses at target depth
# sxx = sigma_xx[idx]
# syy = sigma_yy[idx]
# szz = sigma_zz[idx]
# sxy = sigma_xy[idx]
# Pf_target = Pf[idx]
# xi_target = xi[idx]

# # Principal stresses & orientation
# s_tensor = np.array([[sxx, sxy], [sxy, syy]])
# eigvals, eigvecs = np.linalg.eigh(s_tensor)  # ascending order
# s2, s1 = eigvals  # s1 = major, s2 = minor
# theta_deg = np.degrees(np.arctan2(eigvecs[1, 0], eigvecs[0, 0]))  # angle of s1

# # Report
# print(f"Depth (array snap-to):  {depth_exact / 1e3:.3f} km")
# print(f"Pore-fluid pressure:   {Pf_target / 1e6:8.2f}  MPa")
# print(f"szz (vertical):        {szz / 1e6:8.2f}  MPa")
# print(
#     f"sxx, syy, sxy:         {sxx / 1e6:8.2f}, {syy / 1e6:8.2f}, {sxy / 1e6:8.2f} MPa"
# )
# print(f"s1 (major):            {s1 / 1e6:8.2f}  MPa")
# print(f"s2 (minor):            {s2 / 1e6:8.2f}  MPa")
# print(f"Orientation of s1:     {theta_deg:6.2f}°  (CCW from +x)")
# print(f"Invariant xi:          {xi_target:.3f}")

# ------------------------
# Middle subplot: Effective Stress Components
# ------------------------
# plt.plot(Pf / 1e6, depths / 1e3, label=r"$P_f$")
ax2.plot(np.abs(sigma_zz) / 1e6, depths / 1e3, label=r"$\sigma_{zz}$")
ax2.plot(np.abs(sigma_xx) / 1e6, depths / 1e3, label=r"$\sigma_{xx}$")
ax2.plot(np.abs(sigma_yy) / 1e6, depths / 1e3, label=r"$\sigma_{yy}$")
ax2.plot(sigma_xy / 1e6, depths / 1e3, label=r"$\sigma_{xy}$")
ax2.plot(
    static_shear_strength / 1e6,
    depths / 1e3,
    label="Static Shear Strength",
    linestyle="--",
    alpha=0.8,
)
ax2.plot(
    residual_shear_strength / 1e6,
    depths / 1e3,
    label="Residual Shear Strength",
    linestyle="--",
    alpha=0.8,
)

ax2.invert_yaxis()
ax2.set_ylabel("Depth (km)", fontsize=20)
ax2.set_xlabel("Stress (MPa)", fontsize=20)
ax2.set_title("Effective Stress Components vs Depth", fontsize=20)
ax2.legend(loc="best", fontsize=10)
ax2.grid(True, which="both", linestyle=":")
for lab in ax2.get_xticklabels() + ax2.get_yticklabels():
    lab.set_fontproperties(tick_prop)

# Third subplot: Strain Invariants
# ax3.plot(I1, depths / 1e3, label=r"$I_1$")
# ax3.plot(I2, depths / 1e3, label=r"$I_2$")
ax3.plot(xi, depths / 1e3, label=r"$\xi$")
ax3.invert_yaxis()
ax3.set_ylabel("Depth (km)", fontsize=20)
ax3.set_xlabel("Invariant values", fontsize=20)
ax3.set_title("Strain Invariants vs Depth", fontsize=20)
ax3.set_xlim(-1.8, 0)
ax3.legend(loc="best", fontsize=10)
ax3.grid(True, which="both", linestyle=":")
for lab in ax3.get_xticklabels() + ax3.get_yticklabels():
    lab.set_fontproperties(tick_prop)

# Save combined figure
plt.tight_layout()
plt.savefig("stress_and_strain_combined.png", dpi=300)
plt.show()
