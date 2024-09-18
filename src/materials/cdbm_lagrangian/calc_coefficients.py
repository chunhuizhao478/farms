from sympy import symbols, Eq, solve

# Define the unknowns
a0, a1, a2, a3 = symbols('a0 a1 a2 a3')

# Define the known quantities
xi1, xi_d, chi, mu, gamma_r, lambda_o, mu_o, xi_o = symbols('xi1 xi_d chi mu gamma_r lambda_o mu_o xi_o')

# Define the equations
eq1 = Eq(2 * a2 + a1 / xi1 + 3 * a3 * xi1, 0)
eq2 = Eq(2 * a0 + a1 * xi1 - a3 * xi1**3, 0)
eq3 = Eq(a0, chi * mu)
eq4 = Eq(a0 + a1 * xi_d + a2 * xi_d**2 + a3 * xi_d**3, (mu_o + xi_o * gamma_r) - gamma_r * xi_d + (lambda_o / 2) * xi_d**2)

# Solve the system of equations for a0, a1, a2, and a3
solution = solve([eq1, eq2, eq3, eq4], (a0, a1, a2, a3))

#####