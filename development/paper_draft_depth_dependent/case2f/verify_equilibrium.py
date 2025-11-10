"""
Verify equilibrium between InitialStressStrainTPV26 and EffectiveBodyForceTPV26
Check that dσ'_zz/dz + f_z = 0 for all regions
"""

import numpy as np
import matplotlib.pyplot as plt

# Parameters from input file
fluid_density = 1000  # kg/m³
rock_density = 2670   # kg/m³
gravity = 9.8         # m/s²
lambda_pp = 0.9       # pore pressure ratio
A = 6000              # transition start depth (m)
B = 8000              # transition end depth (m)

# Create depth array
depths = np.linspace(0, 20000, 1000)

# Compute pore pressure and its gradient
Pf = np.zeros_like(depths)
dPf_dz = np.zeros_like(depths)

for i, z in enumerate(depths):
    if z <= A:
        # Region 1: Hydrostatic
        Pf[i] = fluid_density * gravity * z
        dPf_dz[i] = fluid_density * gravity
    elif z > A and z <= B:
        # Region 2: Quadratic transition
        Pf_A = fluid_density * gravity * A
        Pf_B_target = lambda_pp * rock_density * gravity * B
        s = (z - A) / (B - A)
        Pf[i] = Pf_A + (Pf_B_target - Pf_A) * s**2
        # dPf/dz = 2 * (Pf_B_target - Pf_A) / (B - A) * s
        dPf_dz[i] = 2.0 * (Pf_B_target - Pf_A) / (B - A) * s
    else:
        # Region 3: Scaled lithostatic
        Pf[i] = lambda_pp * rock_density * gravity * z
        dPf_dz[i] = lambda_pp * rock_density * gravity

# Compute vertical stress and its gradient
sigma_zz = -rock_density * gravity * depths
d_sigma_zz_dz = -rock_density * gravity * np.ones_like(depths)

# Compute effective vertical stress and its gradient
sigma_eff_zz = sigma_zz + Pf
d_sigma_eff_zz_dz = d_sigma_zz_dz + dPf_dz

# Compute body force from EffectiveBodyForceTPV26
# In the kernel: f_eff = -1.0 * (rho*g - dPf/dz), then residual = -f_eff * test
# So the effective body force applied is: (rho*g - dPf/dz) (positive downward)
f_z_downward = rock_density * gravity - dPf_dz  # positive downward

# Check equilibrium: dσ'_zz/dz + f_z = 0 (both use same sign convention)
# Since σ'_zz is compression negative, dσ'_zz/dz is negative downward
# So equilibrium is: dσ'_zz/dz + f_z_downward = 0
equilibrium_residual = d_sigma_eff_zz_dz + f_z_downward

# Plot results
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Plot 1: Pore pressure and gradient
ax1 = axes[0, 0]
ax1_twin = ax1.twinx()
ax1.plot(depths/1000, Pf/1e6, 'b-', label='Pf', linewidth=2)
ax1_twin.plot(depths/1000, dPf_dz/1e3, 'r--', label='dPf/dz', linewidth=2)
ax1.axvline(A/1000, color='k', linestyle=':', alpha=0.5, label='A=6km')
ax1.axvline(B/1000, color='k', linestyle=':', alpha=0.5, label='B=8km')
ax1.set_xlabel('Depth (km)', fontsize=12)
ax1.set_ylabel('Pore Pressure (MPa)', fontsize=12, color='b')
ax1_twin.set_ylabel('dPf/dz (kPa/m)', fontsize=12, color='r')
ax1.tick_params(axis='y', labelcolor='b')
ax1_twin.tick_params(axis='y', labelcolor='r')
ax1.set_title('Pore Pressure and Gradient', fontsize=14)
ax1.grid(True, alpha=0.3)
ax1.legend(loc='upper left')

# Plot 2: Effective stress and gradient
ax2 = axes[0, 1]
ax2_twin = ax2.twinx()
ax2.plot(depths/1000, -sigma_eff_zz/1e6, 'b-', label="σ'_zz", linewidth=2)
ax2_twin.plot(depths/1000, -d_sigma_eff_zz_dz/1e3, 'r--', label="dσ'_zz/dz", linewidth=2)
ax2.axvline(A/1000, color='k', linestyle=':', alpha=0.5)
ax2.axvline(B/1000, color='k', linestyle=':', alpha=0.5)
ax2.set_xlabel('Depth (km)', fontsize=12)
ax2.set_ylabel("Effective Stress σ'_zz (MPa)", fontsize=12, color='b')
ax2_twin.set_ylabel("dσ'_zz/dz (kPa/m)", fontsize=12, color='r')
ax2.tick_params(axis='y', labelcolor='b')
ax2_twin.tick_params(axis='y', labelcolor='r')
ax2.set_title('Effective Vertical Stress and Gradient', fontsize=14)
ax2.grid(True, alpha=0.3)
ax2.legend(loc='upper left')

# Plot 3: Body force
ax3 = axes[1, 0]
ax3.plot(depths/1000, f_z_downward/1e3, 'g-', linewidth=2, label='Body Force f_z (downward)')
ax3.axvline(A/1000, color='k', linestyle=':', alpha=0.5)
ax3.axvline(B/1000, color='k', linestyle=':', alpha=0.5)
ax3.axhline(0, color='k', linestyle='-', alpha=0.3)
ax3.set_xlabel('Depth (km)', fontsize=12)
ax3.set_ylabel('Body Force f_z (kPa/m³)', fontsize=12)
ax3.set_title('Effective Body Force (positive = downward)', fontsize=14)
ax3.grid(True, alpha=0.3)
ax3.legend()

# Plot 4: Equilibrium residual (should be zero everywhere)
ax4 = axes[1, 1]
ax4.plot(depths/1000, equilibrium_residual, 'r-', linewidth=2, label="dσ'_zz/dz + f_z")
ax4.axvline(A/1000, color='k', linestyle=':', alpha=0.5)
ax4.axvline(B/1000, color='k', linestyle=':', alpha=0.5)
ax4.axhline(0, color='g', linestyle='--', linewidth=2, alpha=0.7, label='Perfect Equilibrium')
ax4.set_xlabel('Depth (km)', fontsize=12)
ax4.set_ylabel('Equilibrium Residual (Pa/m)', fontsize=12)
ax4.set_title('Equilibrium Check: dσ\'_zz/dz + f_z', fontsize=14)
ax4.grid(True, alpha=0.3)
ax4.legend()
ax4.set_ylim([-1, 1])  # Zoom in to see any non-zero values

plt.tight_layout()
plt.savefig('equilibrium_verification.png', dpi=300)
plt.show()

# Print statistics
print("\n" + "="*70)
print("EQUILIBRIUM VERIFICATION FOR lambda_pp = {:.2f}".format(lambda_pp))
print("="*70)
print("\nRegion 1 (0-6 km): Hydrostatic")
mask1 = depths <= A
print(f"  Effective body force: {f_z_downward[mask1][0]/1e3:.2f} kPa/m³ (constant)")
print(f"  dσ'_zz/dz: {d_sigma_eff_zz_dz[mask1][0]/1e3:.2f} kPa/m³")
print(f"  Equilibrium residual: max = {np.max(np.abs(equilibrium_residual[mask1])):.2e} Pa/m")

print("\nRegion 2 (6-8 km): Quadratic Transition")
mask2 = (depths > A) & (depths <= B)
print(f"  Effective body force range: {f_z_downward[mask2].min()/1e3:.2f} to {f_z_downward[mask2].max()/1e3:.2f} kPa/m³")
print(f"  dσ'_zz/dz range: {d_sigma_eff_zz_dz[mask2].min()/1e3:.2f} to {d_sigma_eff_zz_dz[mask2].max()/1e3:.2f} kPa/m³")
print(f"  Equilibrium residual: max = {np.max(np.abs(equilibrium_residual[mask2])):.2e} Pa/m")

print("\nRegion 3 (>8 km): Scaled Lithostatic (lambda_pp = {:.2f})".format(lambda_pp))
mask3 = depths > B
print(f"  Effective body force: {f_z_downward[mask3][0]/1e3:.2f} kPa/m³ (constant)")
print(f"  dσ'_zz/dz: {d_sigma_eff_zz_dz[mask3][0]/1e3:.2f} kPa/m³")
print(f"  Effective stress fraction: {(1-lambda_pp)*100:.1f}% of overburden")
print(f"  Equilibrium residual: max = {np.max(np.abs(equilibrium_residual[mask3])):.2e} Pa/m")

print("\nGLOBAL EQUILIBRIUM:")
print(f"  Maximum equilibrium residual: {np.max(np.abs(equilibrium_residual)):.2e} Pa/m")
print(f"  Mean equilibrium residual: {np.mean(np.abs(equilibrium_residual)):.2e} Pa/m")

if np.max(np.abs(equilibrium_residual)) < 1e-10:
    print("\n✓ EQUILIBRIUM IS SATISFIED (residual < 1e-10 Pa/m)")
else:
    print("\n✗ WARNING: Equilibrium residual is not zero!")

print("="*70 + "\n")
