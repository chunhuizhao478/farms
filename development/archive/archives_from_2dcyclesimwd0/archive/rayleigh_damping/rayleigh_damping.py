import numpy as np

# Given data
zeta = 0.1 # desired damping ratio
fmin = 19  # Hz
fmax = 39  # Hz

omega_min = 2 * np.pi * fmin  # rad/s
omega_max = 2 * np.pi * fmax  # rad/s

# Calculation
A = np.array([[1/omega_min, omega_min],
              [1/omega_max, omega_max]])
b = np.array([2*zeta, 2*zeta])
alpha, beta = np.linalg.solve(A, b)

print("alpha : ", alpha, "beta : ", beta)