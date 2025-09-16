import numpy as np

# Given parameters
xi_o = -0.8
lmbda_mu_ratio = 1  # lambda_o / mu_o = 1 based on common assumption unless otherwise provided

# Solve for q from equation (A.7)
# xi_o = -sqrt(3) / sqrt(2 * q^2 * (lmbda_mu_ratio + 2/3)^2 + 1)

# Rearranging:
# xi_o^2 = 3 / (2 * q^2 * (lmbda_mu_ratio + 2/3)^2 + 1)
# (2 * q^2 * (lmbda_mu_ratio + 2/3)^2 + 1) = 3 / xi_o^2
# 2 * q^2 * (lmbda_mu_ratio + 2/3)^2 = 3 / xi_o^2 - 1
# q^2 = (3 / xi_o^2 - 1) / (2 * (lmbda_mu_ratio + 2/3)^2)

numerator = 3 / xi_o**2 - 1
denominator = 2 * (lmbda_mu_ratio + 2/3)**2
q_squared = numerator / denominator
q = np.sqrt(q_squared)

# Solve for phi from equation (A.8)
# q = sin(phi) / (1 - sin(phi)/3)
# Rearranging: q * (1 - sin(phi)/3) = sin(phi)
# q - q * sin(phi)/3 = sin(phi)
# q = sin(phi) * (1 + q/3)
# sin(phi) = q / (1 + q/3)

sin_phi = q / (1 + q / 3)
phi_rad = np.arcsin(sin_phi)
phi_deg = np.degrees(phi_rad)

print(f"Angle of internal friction (phi): {phi_deg:.2f} degrees")