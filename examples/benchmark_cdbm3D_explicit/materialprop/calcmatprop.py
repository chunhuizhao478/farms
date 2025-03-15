# Reinitialize values after execution state reset
import numpy as np

# Given new values
mu = 32.04e9   # Shear modulus (Pa)
lambda_ = 32.04e9  # Lame's first parameter (Pa)
rho = 2700  # Density (kg/m^3)

# Compute shear wave speed Vs
Vs = np.sqrt(mu / rho)

# Compute Young's modulus using given Lame parameters
E = mu * (3 * lambda_ + 2 * mu) / (lambda_ + mu)

# Compute Poisson's ratio
nu = lambda_ / (2 * (lambda_ + mu))

# Compute first compressional wave speed formula Vp (including Poisson's ratio effect)
Vp_full = np.sqrt((E * (1 - nu)) / (rho * (1 + nu) * (1 - 2 * nu)))

# Compute second compressional wave speed formula Vp (using Lame parameters)
Vp_lame = np.sqrt((lambda_ + 2 * mu) / rho)

# Display results
print(f"Young's modulus E: {E:.2f} Pa")
print(f"Poisson's ratio nu: {nu:.2f}")
print(f'Shear wave speed Vs: {Vs:.2f} m/s')
print(f'First compressional wave speed Vp (full formula): {Vp_full:.2f} m/s')
print(f'Second compressional wave speed Vp (using Lame parameters): {Vp_lame:.2f} m/s')
