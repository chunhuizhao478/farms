import numpy as np
import matplotlib.pyplot as plt
from sympy import symbols, Eq, solve

# Define the function for the negative root solution
def alpha_func_root1(xi, lambda_val=1, mu_o_val=1, xi_o_val=0.5, gamma_r_val=1):
    term1 = lambda_val * xi**3 - 6 * lambda_val * xi_o_val + 6 * mu_o_val * xi - 8 * mu_o_val * xi_o_val
    term2 = np.sqrt(lambda_val**2 * xi**6 - 12 * lambda_val**2 * xi**3 * xi_o_val + 36 * lambda_val**2 * xi_o_val**2
                    + 12 * lambda_val * mu_o_val * xi**4 - 16 * lambda_val * mu_o_val * xi**3 * xi_o_val
                    - 72 * lambda_val * mu_o_val * xi**2 + 72 * lambda_val * mu_o_val * xi * xi_o_val
                    + 72 * lambda_val * mu_o_val - 12 * mu_o_val**2 * xi**2 + 48 * mu_o_val**2)
    denominator = 2 * gamma_r_val * (3 * xi**2 - 6 * xi * xi_o_val + 4 * xi_o_val**2 - 3)
    return (term1 - term2) / denominator

lambda_o = 30e9
mu_o = 30e9
xi_o = -0.8
xi_d = -0.9
chi = 0.8
xig = 0

# The expression for gamma_r:
gamma_r = (
    -xi_o * (-xi_o**2 + 6*lambda_o + 2*mu_o)
    - np.sqrt(
        (xi_o**2 + 2*mu_o)
        * (
            lambda_o * (xi_o**2 - 12*xi_o + 36)
            - 6*xi_o*mu_o
            + 24*mu_o
        )
    )
) / (2 * (xi_o**2 - 3))

# For reference, xi1 in your snippet:
xi1 = xi_o + np.sqrt(xi_o**2 + 2 * mu_o / lambda_o)

print("gamma_r =", gamma_r)
print("xi1 =", xi1)

#get alphacr(xi_given)
alpha_g = alpha_func_root1(xig, lambda_o, mu_o, xi_o, gamma_r)
alpha_1 = alpha_func_root1(xi1, lambda_o, mu_o, xi_o, gamma_r)

# Define the unknowns
a0, a1, a2, a3 = symbols('a0 a1 a2 a3')

# Define the known quantities
# xi1, xi_d, chi, mu, gamma_r, lambda_o, mu_o, xi_o = symbols('xi1 xi_d chi mu gamma_r lambda_o mu_o xi_o')

# Define the equations
eq1 = Eq(2 * a2 + a1 / xi1 + 3 * a3 * xi1, 0)
eq2 = Eq(2 * a0 + a1 * xi1 - a3 * xi1**3, 0)
eq3 = Eq(a0 + a1 * xig  + a2 * xig**2 + a3 * xig**3, chi * ((mu_o + alpha_g * xi_o * gamma_r) - alpha_g * gamma_r * xig + (lambda_o / 2) * xig**2))
eq4 = Eq(a0 + a1 * xi_d + a2 * xi_d**2 + a3 * xi_d**3, (mu_o + xi_d * gamma_r) - gamma_r * xi_d + (lambda_o / 2) * xi_d**2)

# eq1 = Eq(a0 - chi * (mu_o + alpha_1 * xi_o * gamma_r), 0)
# eq2 = Eq(a1 + chi * (alpha_1 * gamma_r), 0)
# eq3 = Eq(a2 - chi * (lambda_o / 2), 0)
# eq4 = Eq(a0 + a1 * xi_d + a2 * xi_d**2 + a3 * (xi_d-xi1)**3, (mu_o + xi_d * gamma_r) - gamma_r * xi_d + (lambda_o / 2) * xi_d**2)

# Solve the system of equations for a0, a1, a2, and a3
solution = solve([eq1, eq2, eq3, eq4], (a0, a1, a2, a3))

# Extract the solution
sol = solution
a0 = sol[a0]
a1 = sol[a1]
a2 = sol[a2]
a3 = sol[a3]

# Print the solution
print(solution)

alpha = [0.6,0.7,0.8,0.9,1.0]
mu = np.zeros(len(alpha))
gamma = np.zeros(len(alpha))

for i in range(len(alpha)):
    mu[i] = mu_o + alpha[i] * xi_o * gamma_r
    gamma[i] = gamma_r * alpha[i]

xi = np.linspace(-1.7,1.7,100)

shear_modulus_elastic = np.zeros((len(alpha),len(xi)))
shear_modulus_granular = np.zeros(len(xi))

for i in range(len(alpha)):
    for j in range(len(xi)):
        shear_modulus_elastic[i,j] = (2 * mu[i] - gamma[i] * xi[j])/2
    shear_modulus_granular = (2 * a0 + a1 * xi - a3 * xi**3 )/2
    # shear_modulus_granular = (2 * a0 + a1 * xi - a3)/2

plt.figure()
plt.plot(xi, shear_modulus_elastic[0,:]/1e9, label=r'$\alpha = 0.6$')
plt.plot(xi, shear_modulus_elastic[1,:]/1e9, label=r'$\alpha = 0.7$')
plt.plot(xi, shear_modulus_elastic[2,:]/1e9, label=r'$\alpha = 0.8$')
plt.plot(xi, shear_modulus_elastic[3,:]/1e9, label=r'$\alpha = 0.9$')
plt.plot(xi, shear_modulus_elastic[4,:]/1e9, label=r'$\alpha = 1.0$')
plt.plot(xi, shear_modulus_granular/1e9, 'k-*', label='granular')
plt.xlabel(r'$\xi$')
plt.ylabel(r'$G (MPa)$')
plt.xlim([-1.7,1.7])
# plt.ylim([-1,15])
plt.grid()
plt.legend()
plt.show()