import numpy as np
import matplotlib.pyplot as plt

# Given stress components
sigma_xx = -50e6  # Pa
sigma_yy = -50e6  # Pa
tau_xy = 13e6     # Pa

print(f"Input stresses:")
print(f"σ_xx = {sigma_xx/1e6:.1f} MPa")
print(f"σ_yy = {sigma_yy/1e6:.1f} MPa")
print(f"τ_xy = {tau_xy/1e6:.1f} MPa")

# Compute Mohr circle center and radius
C = (sigma_xx + sigma_yy) / 2  # Center on σ-axis
R = np.sqrt(((sigma_xx - sigma_yy) / 2)**2 + tau_xy**2)  # Radius

print(f"\nMohr Circle Parameters:")
print(f"Center: {C/1e6:.1f} MPa")
print(f"Radius: {R/1e6:.1f} MPa")

# Calculate principal stresses
sigma_1 = C + R  # Maximum principal stress
sigma_2 = C - R  # Minimum principal stress

print(f"\nPrincipal Stresses:")
print(f"σ₁ (max) = {sigma_1/1e6:.1f} MPa")
print(f"σ₂ (min) = {sigma_2/1e6:.1f} MPa")

# Calculate maximum shear stress
tau_max = R

print(f"τ_max = {tau_max/1e6:.1f} MPa")

# Parametric angles for circle
theta = np.linspace(0, 2 * np.pi, 400)
sigma_circle = C + R * np.cos(theta)
tau_circle = R * np.sin(theta)

# Failure envelope parameters
phi = np.deg2rad(45)  # friction angle
c = 0  # cohesion
sigma_n = np.linspace(-80e6, 20e6, 400)

# Mohr-Coulomb failure envelope: τ = c + σ_n * tan(φ)
# For compression (negative normal stress), we need the absolute value
tau_env_upper = c + np.abs(sigma_n) * np.tan(phi)
tau_env_lower = -(c + np.abs(sigma_n) * np.tan(phi))

# Plotting
plt.figure(figsize=(10, 8))

# Plot Mohr circle
plt.plot(sigma_circle/1e6, tau_circle/1e6, 'b-', linewidth=2, label="Mohr Circle")

# Plot failure envelope (both upper and lower)
plt.plot(sigma_n/1e6, tau_env_upper/1e6, 'r-', linewidth=2, label="Failure Envelope (φ=45°)")
plt.plot(sigma_n/1e6, tau_env_lower/1e6, 'r-', linewidth=2)

# Plot stress state points
plt.scatter([sigma_xx/1e6, sigma_yy/1e6], [tau_xy/1e6, -tau_xy/1e6], 
           color='green', s=100, marker='o', label="Stress State Points", zorder=5)

# Plot principal stress points
plt.scatter([sigma_1/1e6, sigma_2/1e6], [0, 0], 
           color='red', s=100, marker='x', label="Principal Stresses", zorder=5)

# Plot center
plt.scatter([C/1e6], [0], color='blue', s=50, marker='+', label="Circle Center", zorder=5)

# Add grid and axes
plt.axhline(0, color='k', linewidth=0.5, alpha=0.5)
plt.axvline(0, color='k', linewidth=0.5, alpha=0.5)
plt.grid(True, alpha=0.3)

# Labels and formatting
plt.xlabel('Normal Stress σ (MPa)', fontsize=12)
plt.ylabel('Shear Stress τ (MPa)', fontsize=12)
plt.title('Mohr Circle Analysis with Mohr-Coulomb Failure Envelope\n(φ=45°, c=0)', fontsize=14)
plt.legend(loc='upper right')

# Equal aspect ratio for proper circle appearance
plt.axis('equal')
plt.tight_layout()

# Check if failure occurs
circle_intersects_envelope = False
for i in range(len(sigma_circle)):
    sigma_val = sigma_circle[i]
    tau_val = abs(tau_circle[i])
    tau_envelope_val = c + abs(sigma_val) * np.tan(phi)
    if tau_val >= tau_envelope_val:
        circle_intersects_envelope = True
        break

print(f"\nFailure Analysis:")
if circle_intersects_envelope:
    print("FAILURE PREDICTED - Mohr circle intersects failure envelope")
else:
    print("NO FAILURE - Mohr circle is within failure envelope")

plt.show()