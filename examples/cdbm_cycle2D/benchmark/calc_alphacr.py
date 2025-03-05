from sympy import symbols, Eq, solve

# Define the symbols
alpha, xi, gamma_r, mu_o, xi_o, lambda_ = symbols('alpha xi gamma_r mu_o xi_o lambda_')

#
mu = mu_o + alpha * xi_o * gamma_r
gamma = alpha * gamma_r

# Define the expanded quadratic equation terms (simplified)
eq1 = ( 2 * mu - gamma * xi ) ** 2 + ( 2 * mu - gamma * xi ) * ( 3 * lambda_ - gamma * xi ) + ( lambda_ * gamma * xi - gamma ** 2 ) * ( 3 - gamma ** 2 )
eq2 = ( 2 * mu - gamma * xi )

# Combine all terms and set equal to zero
equation1 = Eq(eq1, 0)
equation2 = Eq(eq2, 0)

# Solve for alpha
alpha_solution1 = solve(equation1, alpha)
alpha_solution2 = solve(equation2, alpha)

#print(alpha_solution1)
#print(alpha_solution2)
####

# Full code for computing and plotting the negative root solution

import numpy as np
import matplotlib.pyplot as plt

# Define the function for the negative root solution
def alpha_func_root1(xi, lambda_val=1, mu_o_val=1, xi_o_val=0.5, gamma_r_val=1):
    term1 = lambda_val * xi**3 - 6 * lambda_val * xi_o_val + 6 * mu_o_val * xi - 8 * mu_o_val * xi_o_val
    term2 = np.sqrt(lambda_val**2 * xi**6 - 12 * lambda_val**2 * xi**3 * xi_o_val + 36 * lambda_val**2 * xi_o_val**2
                    + 12 * lambda_val * mu_o_val * xi**4 - 16 * lambda_val * mu_o_val * xi**3 * xi_o_val
                    - 72 * lambda_val * mu_o_val * xi**2 + 72 * lambda_val * mu_o_val * xi * xi_o_val
                    + 72 * lambda_val * mu_o_val - 12 * mu_o_val**2 * xi**2 + 48 * mu_o_val**2)
    denominator = 2 * gamma_r_val * (3 * xi**2 - 6 * xi * xi_o_val + 4 * xi_o_val**2 - 3)
    return (term1 - term2) / denominator

def alpha_func_root2(xi, mu_o_val=1, xi_o_val=0.5, gamma_r_val=1):
    return 2*mu_o_val/(gamma_r_val*(xi - 2*xi_o_val))

# Define the given parameter values
lambda_val = 30e9
mu_o_val = 30e9
xi_o_val = -0.8
gamma_r_val = 34.785e9

# Generate xi values from -sqrt(3) to sqrt(3)
xi_vals_plot = np.linspace(-np.sqrt(3), np.sqrt(3), 400)

# Compute the alpha values for the negative root with the given parameters
alpha_vals_root_1 = alpha_func_root1(xi_vals_plot, lambda_val, mu_o_val, xi_o_val, gamma_r_val)
alpha_vals_root_2 = alpha_func_root2(xi_vals_plot, mu_o_val, xi_o_val, gamma_r_val)

#print xi = 0
print(alpha_func_root1(0, lambda_val, mu_o_val, xi_o_val, gamma_r_val))
exit()

# Restrict alpha between 0 and 1
alpha_vals_root_1_restricted = np.clip(alpha_vals_root_1, 0, 1)
alpha_vals_root_2_restricted = np.clip(alpha_vals_root_2, 0, 1)

# Line plot for the negative root solution with custom parameters
plt.figure(figsize=(8, 6))
plt.plot(xi_vals_plot, alpha_vals_root_1_restricted, label=r'$\alpha$ (negative root, restricted to [0, 1])', color='r')
plt.plot(xi_vals_plot, alpha_vals_root_2_restricted, label=r'$\alpha$ (negative root, restricted to [0, 1])', color='b')
plt.title(r'Negative root solution for $\alpha$ as a function of $\xi$ (restricted to [0, 1], custom parameters)')
plt.xlabel(r'$\xi$')
plt.ylabel(r'$\alpha$')
plt.grid(True)
plt.legend()
plt.show()