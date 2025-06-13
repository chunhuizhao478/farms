#!/usr/bin/env python3

import sympy

def main():
    # --------------------------------------------------
    # 1) Define symbols
    # --------------------------------------------------
    alpha = sympy.Symbol('alpha', real=True)  # a0 = alpha (free parameter)

    # Unknowns in the system (we will solve for a1, a2, a3 in terms of alpha)
    a1, a2, a3 = sympy.symbols('a1 a2 a3', real=True)

    # Constants from the original system of 3 equations:
    #   xi1, xi_d, mu0, xi0, gamma, r, lam
    # plus extra 'mu' and 'xi' for the inequality
    xi1, xi_d = sympy.symbols('xi1 xi_d', positive=True)  # assume positive
    mu0, xi0, gamma_r, lam = sympy.symbols('mu0 xi0 gamma_r lam', real=True)
    mu, xi = sympy.symbols('mu xi', real=True)

    # --------------------------------------------------
    # 2) Write down the 3 equations
    #    (I)   2 a2 + a1/xi1 + 3 a3 xi1 = 0
    #    (II)  2 a0 + a1 xi1 - a3 xi1^3 = 0
    #    (III) a0 + a1 xi_d + a2 xi_d^2 + a3 xi_d^3 = R
    # where a0 = alpha (the free parameter).
    # --------------------------------------------------

    # Define R for convenience:
    # R = mu0 + xi0*gamma - gamma*r*xi_d + (lam/2)*xi_d^2
    R = mu0 + xi0*gamma_r - gamma_r*xi_d + (lam/sympy.Integer(2))*xi_d**2

    eq1 = sympy.Eq(2*a2 + a1/xi1 + 3*a3*xi1, 0)
    eq2 = sympy.Eq(2*alpha + a1*xi1 - a3*xi1**3, 0)
    eq3 = sympy.Eq(alpha + a1*xi_d + a2*xi_d**2 + a3*xi_d**3, R)

    # --------------------------------------------------
    # 3) Solve for a1, a2, a3 in terms of alpha
    # --------------------------------------------------
    sol = sympy.solve([eq1, eq2, eq3], [a1, a2, a3], dict=True)
    # Typically we get one solution in 'sol'
    a1_expr = sol[0][a1]
    a2_expr = sol[0][a2]
    a3_expr = sol[0][a3]

    # --------------------------------------------------
    # 4) Define the inequality f(alpha) < 0:
    #
    # f(alpha) = (
    #     (alpha - mu)
    #     + (a1_expr + gamma_r)*xi
    #     + (a2_expr - lam_over_2)*xi**2
    #     + a3_expr*xi**3
    #
    # --------------------------------------------------
    lam_over_2 = lam / 2
    f_expr = (
        (alpha - (mu0 - xi0*gamma_r))
        + (a1_expr + gamma_r)*xi
        + (a2_expr - lam_over_2)*xi**2
        + a3_expr*xi**3
    )

    # --------------------------------------------------
    # 5) Substitute numerical values and solve f_expr(alpha) < 0
    # --------------------------------------------------
    # Example numeric assignments (change these as needed):
    numeric_subs = {
        xi1: 0.8248,
        xi_d: -0.9,
        mu0: 30e9,
        xi0: -0.8,
        gamma_r : 3.4785e10,
        lam: 30e9,   # <-- 'lam' for both R and the inequality
        mu: 30e9,    # 'mu' in the inequality
        xi: 0,    # the xi in the inequality expression
    }

    # Create a purely symbolic expression in alpha by substituting numeric values
    f_num = f_expr.subs(numeric_subs)

    # We'll solve the inequality f_num(alpha) < 0.
    # In Sympy, we can use 'solveset' with StrictLessThan:
    alpha_solution = sympy.solveset(sympy.StrictLessThan(f_num, 0), alpha, sympy.Reals)

    # --------------------------------------------------
    # 6) Print results
    # --------------------------------------------------
    print("\n=== Parametric Solutions for a1, a2, a3 in Terms of alpha ===")
    print("a1(alpha) =", a1_expr)
    print("a2(alpha) =", a2_expr)
    print("a3(alpha) =", a3_expr)

    print("\n=== The inequality f(alpha) < 0 ===")
    print("f(alpha) =", f_expr)

    print("\n=== Numeric intervals for alpha that satisfy f_num(alpha) < 0 ===")
    print("Substituting numeric parameters:", numeric_subs)
    print("Solution for alpha in real intervals =")
    print(alpha_solution)

if __name__ == "__main__":
    main()