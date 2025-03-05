from sympy import symbols, Eq, solve

# Define symbols
mu_0, gamma_r, lambda_, xi_o, alpha = symbols('mu_0 gamma_r lambda_ xi_o alpha')

# Assume alpha = 1 and xi = xi_o
alpha_val = 1
xi_val = xi_o

# Define the equation based on the image and equate to 0
equation = Eq((2 * (mu_0 + 1 * xi_o * gamma_r) - gamma_r * xi_o)**2 + (2 * (mu_0 + 1 * xi_o * gamma_r) - (gamma_r) * xi_o) * (3 * lambda_ - gamma_r * xi_o) + (lambda_ * gamma_r * xi_o - gamma_r**2) * (3 - xi_o**2), 0)

# Solve for gamma
gamma_solution = solve(equation, gamma_r)

print(gamma_solution)

####
# Define numerical values for the parameters
xi_o_val = -0.8
lambda_val = 30e9
mu_o_val = 30e9

# Define gamma_r expressions
# use gamma_r1
gamma_r1 = (-xi_o_val * (-lambda_val * xi_o_val**2 + 6 * lambda_val + 2 * mu_o_val) - 
            ((lambda_val * xi_o_val**2 + 2 * mu_o_val) * 
            (lambda_val * xi_o_val**4 - 12 * lambda_val * xi_o_val**2 + 36 * lambda_val - 6 * mu_o_val * xi_o_val**2 + 24 * mu_o_val))**0.5) / (2 * (xi_o_val**2 - 3))

gamma_r2 = (-xi_o_val * (-lambda_val * xi_o_val**2 + 6 * lambda_val + 2 * mu_o_val) + 
            ((lambda_val * xi_o_val**2 + 2 * mu_o_val) * 
            (lambda_val * xi_o_val**4 - 12 * lambda_val * xi_o_val**2 + 36 * lambda_val - 6 * mu_o_val * xi_o_val**2 + 24 * mu_o_val))**0.5) / (2 * (xi_o_val**2 - 3))

print(gamma_r1, gamma_r2)