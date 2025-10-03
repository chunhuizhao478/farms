import math

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


df = pd.read_csv('data_p_tau.csv')

sigma_n = df.iloc[:, 0].to_numpy(dtype=float) / 1e6  # normal stress in MPa
tau = df.iloc[:, 1].to_numpy(dtype=float) / 1e6      # shear stress in MPa
tau2 = df.iloc[:, 2].to_numpy(dtype=float) / 1e6     # second shear stress in MPa

print(f'Total σ_N values read: {len(sigma_n)}')

# Linear regression for the first shear stress
# ----------------------------------------------------------------- #
mask1 = np.isfinite(sigma_n) & np.isfinite(tau)
sigma_tau = sigma_n[mask1]
tau_clean = tau[mask1]

if sigma_tau.size < 2:
    raise SystemExit('Not enough finite data points for first Mohr-Coulomb regression.')

slope, intercept = np.polyfit(sigma_tau, tau_clean, 1)

print(f'Linear regression (dataset 1): τ = {intercept:.3f} + {slope:.3f} σ_N')

cohesion = intercept
phi_rad = math.atan(slope)
phi_deg = math.degrees(phi_rad)

sigma_line = np.linspace(sigma_tau.min(), sigma_tau.max(), 200)
tau_line = intercept + slope * sigma_line

print('Mohr-Coulomb regression (dataset 1): τ = c + σ_N tanφ')
print(f'  c   = {cohesion:.3f} MPa')
print(f'  φ   = {phi_deg:.2f}°')
print(f'  tanφ = {slope:.4f}')

# Linear regression for the second shear stress
# ----------------------------------------------------------------- #
mask2 = np.isfinite(sigma_n) & np.isfinite(tau2)
sigma_tau2 = sigma_n[mask2]
tau2_clean = tau2[mask2]

has_second_fit = sigma_tau2.size >= 2
if has_second_fit:
    slope2, intercept2 = np.polyfit(sigma_tau2, tau2_clean, 1)

    print(f'Linear regression (dataset 2): τ = {intercept2:.3f} + {slope2:.3f} σ_N')

    cohesion2 = intercept2
    phi_rad2 = math.atan(slope2)
    phi_deg2 = math.degrees(phi_rad2)

    sigma_line2 = np.linspace(sigma_tau2.min(), sigma_tau2.max(), 200)
    tau_line2 = intercept2 + slope2 * sigma_line2

    print('Mohr-Coulomb regression (dataset 2): τ = c + σ_N tanφ')
    print(f'  c   = {cohesion2:.3f} MPa')
    print(f'  φ   = {phi_deg2:.2f}°')
    print(f'  tanφ = {slope2:.4f}')
else:
    print('Dataset 2 does not have enough finite values for Mohr-Coulomb regression.')

plt.figure(figsize=(10, 6))
plt.scatter(sigma_tau, tau_clean, color='red', label='$\\hat{\\dot{\\epsilon}} = 10^{-7} 1/s, C_d = 10 1/s$')
if has_second_fit:
    plt.scatter(sigma_tau2, tau2_clean, color='black', label='$\\hat{\\dot{\\epsilon}} = 10^{-5} 1/s, C_d = 10 1/s$')
plt.plot(sigma_line, tau_line, color='red', label='Mohr-Coulomb fit (dataset 1)')
if has_second_fit:
    plt.plot(sigma_line2, tau_line2, color='black', label='Mohr-Coulomb fit (dataset 2)')
plt.xlabel('Confinement Pressure $P$ (MPa)')
plt.ylabel('Peak Stress $\\tau$ (MPa)')
plt.legend()
plt.tight_layout()
plt.show()
