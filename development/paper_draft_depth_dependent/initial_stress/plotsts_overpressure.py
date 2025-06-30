"""Enhanced pore‑pressure and stress profile with a depth‑dependent transition.

Between depths A and B (m), the vertical pore‑pressure gradient switches smoothly
from hydrostatic (density_fluid · g) to lithostatic/over‑pressure (rho · g).
Below B the over‑pressured gradient is retained.  All other stress components
follow the same empirical relationships defined in the original *plotsts.py*:

    σ_xx = b_xx (σ_zz + P_f) − P_f
    σ_yy = b_yy (σ_zz + P_f) − P_f
    σ_xy = b_xy (σ_zz + P_f)

with coefficients valid down to 15.6 km.  Cohesion, static, and residual
strengths are computed exactly as before.

Edit the constants under “USER‑PARAMETERS” to suit your model.
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt

# -------------------------------------------------------------------
# USER‑PARAMETERS
# -------------------------------------------------------------------
A = 8_000.0      # (m) depth where the transition *begins*
B = 12_000.0     # (m) depth where the transition *ends*
z_max = 15_000.0 # (m) modelling depth
n_pts = 300      # number of depth samples

density_fluid = 1_000.0   # (kg/m³) pore‑fluid density
rho           = 2_670.0   # (kg/m³) bulk rock density (lithostatic)
g             = 9.8       # (m/s²) gravitational acceleration

# Stress coefficients (unchanged)
b_xx = 0.926793
b_yy = 1.073206
b_xy = -0.9

# Friction & cohesion parameters (unchanged)
mu_s = 0.85           # static friction coefficient
mu_d = 0.6            # dynamic / residual friction coefficient
c0   = 0.3e6          # Pa, cohesion at 0–4 km
dc_dz = 0.000675e6    # Pa/m, cohesion gradient up to 4 km

# -------------------------------------------------------------------
# DEPTH GRID
# -------------------------------------------------------------------
depths = np.linspace(0.0, z_max, n_pts)  # positive downward (m)

# -------------------------------------------------------------------
# PORE PRESSURE (piecewise‑analytic)
# -------------------------------------------------------------------
delta_rho = rho - density_fluid  # density contrast

Pf = np.empty_like(depths)

# Region 1: Hydrostatic above A
mask1 = depths <= A
Pf[mask1] = density_fluid * g * depths[mask1]

# Region 2: Linear‑gradient transition A < z ≤ B
mask2 = (depths > A) & (depths <= B)
z2 = depths[mask2]
Pf_A = density_fluid * g * A
Pf[mask2] = Pf_A + g * (
    density_fluid * (z2 - A) +
    0.5 * delta_rho * (z2 - A) ** 2 / (B - A)
)

# Region 3: Over‑pressured below B
mask3 = depths > B
z3 = depths[mask3]
Pf_B = density_fluid * g * B + 0.5 * g * delta_rho * (B - A)
Pf[mask3] = Pf_B + rho * g * (z3 - B)

# -------------------------------------------------------------------
# VERTICAL STRESS (compression negative)
# -------------------------------------------------------------------
sigma_zz = -rho * g * depths

# -------------------------------------------------------------------
# HORIZONTAL & SHEAR STRESSES (same empirical relations)
# -------------------------------------------------------------------
mask_coeff = depths <= 15_600.0  # coefficient domain from original file
sigma_xx = np.where(mask_coeff,
                    b_xx * (sigma_zz + Pf) - Pf,
                    sigma_zz)

sigma_yy = np.where(mask_coeff,
                    b_yy * (sigma_zz + Pf) - Pf,
                    sigma_zz)

sigma_xy = np.where(mask_coeff,
                    b_xy * (sigma_zz + Pf),
                    0.0)

# -------------------------------------------------------------------
# COHESION (piecewise as before)
# -------------------------------------------------------------------
mask_cohesion = depths <= 4_000.0
c = np.where(mask_cohesion,
             c0 + dc_dz * (4_000.0 - depths),
             c0)

# -------------------------------------------------------------------
# SHEAR STRENGTHS
# -------------------------------------------------------------------
static_shear_strength    = c + np.abs(mu_s * (sigma_yy + Pf))
residual_shear_strength  = c + np.abs(mu_d * (sigma_yy + Pf))

# -------------------------------------------------------------------
# PLOT
# -------------------------------------------------------------------
plt.figure(figsize=(6, 8))
plt.plot(Pf        / 1e6, depths / 1e3, label=r'$P_f$')
plt.plot(np.abs(sigma_zz) / 1e6, depths / 1e3, label=r'$\sigma_{zz}$')
plt.plot(np.abs(sigma_xx) / 1e6, depths / 1e3, label=r'$\sigma_{xx}$')
plt.plot(np.abs(sigma_yy) / 1e6, depths / 1e3, label=r'$\sigma_{yy}$')
plt.plot(sigma_xy        / 1e6, depths / 1e3, label=r'$\sigma_{xy}$')
plt.plot(static_shear_strength   / 1e6, depths / 1e3,
         label='Static Shear Strength', linestyle='--', color='orange', alpha=0.8)
plt.plot(residual_shear_strength / 1e6, depths / 1e3,
         label='Residual Shear Strength', linestyle='--', color='red', alpha=0.8)

plt.gca().invert_yaxis()
plt.ylabel('Depth (km)')
plt.xlabel('Stress (MPa)')
plt.title('Stress Components vs Depth')
plt.legend(loc='best')
plt.grid(True, which='both', linestyle=':')
plt.tight_layout()
plt.savefig('stress_components_vs_depth.png', dpi=300)
plt.show()

# Elastic constants
mu = 32.04e9              # shear modulus (Pa)
lmbda = 32.04e9           # Lame's first parameter (Pa)
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
I2 = np.einsum('nij,nij->n', strain, strain)
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