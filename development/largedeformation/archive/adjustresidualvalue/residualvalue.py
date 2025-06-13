#!/usr/bin/env python3

def main():
    # Assigned constant values
    xi1    = 0.8248    # e.g., xi1
    xi_d   = -0.9    # e.g., xi_d
    mu0    = 30e9   # e.g., mu0
    xi0    = -0.8    # e.g., xi0
    gamma_r = 3.4785e10   # e.g., gamma
    lambda0 = 30e9   # e.g., lambda0
    
    # Free parameter a0 = alpha
    alpha = 1000e9 #7.41964e9   # choose a value for the free parameter
    
    # Compute R based on the given constants:
    R = mu0 + xi0 * gamma_r - gamma_r * xi_d + (lambda0 / 2.0) * xi_d**2
    print("Computed R =", R)
    
    # Compute a1 using the derived formula:
    denom_a1 = xi_d - 2 * xi_d**2 / xi1 + xi_d**3 / xi1**2
    numer_a1 = R - alpha * (1 - 3 * xi_d**2 / xi1**2 + 2 * xi_d**3 / xi1**3)
    
    a1 = numer_a1 / denom_a1
    
    # Compute a2 and a3:
    a2 = -(3 * alpha + 2 * a1 * xi1) / xi1**2
    a3 = (2 * alpha + a1 * xi1) / xi1**3
    
    # Print the results
    print("\nSolution for the coefficients:")
    print("a0 =", alpha)
    print("a1 =", a1)
    print("a2 =", a2)
    print("a3 =", a3)

    # Compute the granular phase energy
    xi_given = 0.0
    Fb_I2 = alpha + a1 * xi_given + a2 * xi_given**2 + a3 * xi_given**3
    print("\nGranular phase energy Fb_I2 =", Fb_I2/1e8, "x 10^8 J/m^3")


if __name__ == "__main__":
    main()