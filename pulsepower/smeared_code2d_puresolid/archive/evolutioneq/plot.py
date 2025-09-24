import numpy as np
import matplotlib.pyplot as plt

# Updated parameters
paramA = 0.99
paramB = 5000.0
eps0 = 150e6 / 50e9  # strain

# Strain variable (kappa > eps0)
kappa = np.linspace(eps0 * 1.001, 5 * eps0, 500)

# Compute damage omega using updated formula
omega = 1.0 - eps0 / kappa * (1.0 - paramA) - paramA / np.exp(paramB * (kappa - eps0))

# Plot
plt.figure(figsize=(8, 5))
plt.plot(kappa, omega, label=r'$\omega(\kappa)$')
plt.axvline(eps0, color='gray', linestyle='--', label=r'$\varepsilon_0$')
plt.xlabel(r'$\kappa$ (strain)')
plt.ylabel(r'$\omega$ (damage)')
plt.title('Damage Evolution Law (Updated)')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()