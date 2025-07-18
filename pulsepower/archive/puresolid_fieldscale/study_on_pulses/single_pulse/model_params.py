import numpy as np

E = 50e9
nu = 0.373
ft = 25.5e6
Gc = 57
rho = 2600

lch = E * Gc / ( ft * ft )

#AT1
l = 3/8 * lch

# Shear modulus G
G = E / (2 * (1 + nu))
# Bulk modulus K
K = E / (3 * (1 - 2 * nu))
# P-wave speed
Cp = np.sqrt((K + 4.0/3.0 * G) / rho)
# S-wave speed
Cs = np.sqrt(G / rho)

print(f"lch = {lch:.3e} m")
print(f"l = {l:.3e} m")
print(f"Cp = {Cp:.3e} m/s")
print(f"Cs = {Cs:.3e} m/s")