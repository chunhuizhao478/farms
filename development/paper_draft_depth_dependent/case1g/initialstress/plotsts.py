import numpy as np
import matplotlib.pyplot as plt

# Physical constants
density_fluid = 1000      # kg/m^3
rho = 2670                # kg/m^3
g = 9.8                   # m/s^2

# Elastic constants (depth-varying via TPV32 profile)
# Table knots (depth in meters)
depth_knots = np.array([0, 500, 1000, 1600, 2400, 3600, 5000, 9000, 11000, 15000], dtype=float)
vp_knots    = np.array([2200, 3000, 3600, 4400, 4800, 5250, 5500, 5750, 6100, 6300], dtype=float)
vs_knots    = np.array([1050, 1400, 1950, 2500, 2800, 3100, 3250, 3450, 3600, 3700], dtype=float)
rho_knots   = np.array([2200, 2450, 2550, 2600, 2600, 2620, 2650, 2720, 2750, 2900], dtype=float)

# Depth array (m)
depths = np.linspace(0, 20000, 600)  # from surface to 20 km

# Interpolate seismic properties at all depths (piecewise-linear; clamp at ends)
vp  = np.interp(depths, depth_knots, vp_knots,  left=vp_knots[0],  right=vp_knots[-1])
vs  = np.interp(depths, depth_knots, vs_knots,  left=vs_knots[0],  right=vs_knots[-1])
rho_depth = np.interp(depths, depth_knots, rho_knots, left=rho_knots[0], right=rho_knots[-1])

# Derived elastic constants (arrays)
mu_depth   = rho_depth * vs**2
lmbda_depth = np.maximum(rho_depth * vp**2 - 2.0 * mu_depth, 0.0)

# Pore pressure and vertical stress
# Pf remains linear with constant fluid density
Pf = density_fluid * g * depths
# Correct vertical stress: integrate rock density with depth (overburden)
rho_avg = 0.5 * (rho_depth[1:] + rho_depth[:-1])
dz = np.diff(depths)
overburden = np.concatenate(([0.0], np.cumsum(rho_avg * dz)))  # ∫ rho(z) dz
sigma_zz = -g * overburden

# Coefficients for horizontal and shear stresses
b_xx = 0.4
b_yy = 1.073206
b_xy = -0.8

# Piecewise definitions
# Define tapering coefficient Omega(depth)
Omega = np.ones_like(depths)
Omega[(depths > 15000) & (depths <= 20000)] = (20000 - depths[(depths > 15000) & (depths <= 20000)]) / 5000
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
static_shear_strength = c + abs( mu_s * (sigma_yy + Pf) )
mu_d = 0.6
residual_shear_strength = c + abs( mu_d * (sigma_yy + Pf) )

# ------------------------
# Plot Stress Components
# ------------------------
plt.figure(figsize=(6, 8))
plt.plot(Pf/1e6, depths/1e3, label=r'$P_f$')
plt.plot(abs(sigma_zz)/1e6, depths/1e3, label=r'$\sigma_{zz}$') #when plotting, we use positive values for compression
plt.plot(abs(sigma_xx)/1e6, depths/1e3, label=r'$\sigma_{xx}$')
plt.plot(abs(sigma_yy)/1e6, depths/1e3, label=r'$\sigma_{yy}$')
plt.plot(sigma_xy/1e6, depths/1e3, label=r'$\sigma_{xy}$')
plt.plot(static_shear_strength/1e6, depths/1e3, label='Static Shear Strength', linestyle='--', color='orange')
plt.plot(residual_shear_strength/1e6, depths/1e3, label='Residual Shear Strength', linestyle='--', color='red')
plt.gca().invert_yaxis()
plt.ylabel('Depth (km)')
plt.xlabel('Stress (MPa)')
plt.title('Stress Components vs Depth')
plt.legend(loc='best')
plt.grid(True)
plt.tight_layout()
plt.savefig('stress_components_vs_depth.png', dpi=300)
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

# Inverse Hooke's law with depth-varying μ, λ:
#   epsilon_ij = 1/(2μ) s_ij - λ/(2μ(3λ+2μ)) s_kk δ_ij
trace_coeff = -lmbda_depth / (2 * mu_depth * (3 * lmbda_depth + 2 * mu_depth))
strain = np.zeros_like(stress)
enum = np.trace(stress, axis1=1, axis2=2)  # vector of s_kk
# first term (broadcast 1/(2μ(depth)))
strain = stress / (2.0 * mu_depth[:, None, None])
# add diagonal second term with broadcasting
for k in range(3):
    strain[:, k, k] += trace_coeff * enum

# Compute invariants: I1, I2, xi
I1 = np.trace(strain, axis1=1, axis2=2)
I2 = np.einsum('nij,nij->n', strain, strain)
# Stabilize extremely small values near the surface to avoid xi spikes
I2 = np.maximum(I2, 1e-18)
xi = I1 / np.sqrt(I2)

# ------------------------
# Plot Strain Invariants
# ------------------------
plt.figure(figsize=(6, 8))
plt.plot(I1, depths, label=r'$I_1$')
plt.plot(I2, depths, label=r'$I_2$')
plt.plot(xi, depths, label=r'$\xi$')
plt.gca().invert_yaxis()
plt.ylabel('Depth (m)')
plt.xlabel('Invariant values')
plt.title('Strain Invariants vs Depth')
plt.legend(loc='best')
plt.grid(True)
plt.tight_layout()
plt.savefig('strain_invariants_vs_depth.png', dpi=300)
plt.show()

# -----------------------------------------------------------
# Principal stresses and maximum-principal orientation (2-D)
# -----------------------------------------------------------
# --- pick the target depth (m) ---------------------------------------------
depth_target = 7500          # 7.5 km

# ---------------------------------------------------------------------------
# Find the row in `depths` that is closest to the target
idx = int(np.argmin(np.abs(depths - depth_target)))
depth_exact = depths[idx]     # the exact depth value in your array

# ----- stresses already computed in your script ----------------------------
sxx = sigma_xx[idx]
syy = sigma_yy[idx]
szz = sigma_zz[idx]
sxy = sigma_xy[idx]
Pf  = Pf[idx]
xi = xi[idx]  # invariant xi at the target depth

# ----- principal stresses & orientation ------------------------------------
s_tensor = np.array([[sxx, sxy],
                     [sxy, syy]])

eigvals, eigvecs = np.linalg.eigh(s_tensor)   # ascending order
s2, s1 = eigvals                              # s1 = major, s2 = minor
θ_deg = np.degrees(np.arctan2(eigvecs[1, 0],   # angle of s1
                              eigvecs[0, 0]))

# ----- report ---------------------------------------------------------------
print(f"Depth (array snap-to):  {depth_exact/1e3:.3f} km")
print(f"Pore-fluid pressure:   {Pf/1e6:8.2f}  MPa")
print(f"szz (vertical):        {szz/1e6:8.2f}  MPa")
print(f"sxx, syy, sxy:         {sxx/1e6:8.2f}, {syy/1e6:8.2f}, {sxy/1e6:8.2f} MPa")
print(f"s1 (major):            {s1/1e6:8.2f}  MPa")
print(f"s2 (minor):            {s2/1e6:8.2f}  MPa")
print(f"Orientation of s1:     {θ_deg:6.2f}°  (CCW from +x)")
print(f"Invariant xi:          {xi:.3f}")
# ---------------------------------------------------------------------------