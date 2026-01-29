//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>

/**
 * Brent's method for finding roots of a function in a bracketed interval.
 *
 * This is a port of the Fortran ZEROIN algorithm from Netlib, which combines
 * bisection, linear interpolation, and inverse quadratic interpolation for
 * robust and fast root finding.
 *
 * Reference: Brent, R. P. (1973). Algorithms for Minimization without Derivatives.
 * Implementation based on Tandem SEAS code.
 */
class BrentRootFinder
{
public:
  /**
   * Find a root of the function F in the interval [a, b].
   *
   * @param a Left endpoint of the bracketing interval
   * @param b Right endpoint of the bracketing interval
   * @param F Function to find root of (must satisfy F(a) and F(b) have opposite signs)
   * @param tol Tolerance for convergence (default: 0 uses machine epsilon)
   * @return The root x such that F(x) ≈ 0
   * @throws std::runtime_error if F(a) and F(b) have the same sign
   */
  static double zeroIn(double a, double b, std::function<double(double)> F, double tol = 0.0)
  {
    const double eps = std::numeric_limits<double>::epsilon();

    double Fa = F(a);
    double Fb = F(b);

    // Check that we have a valid bracket
    if (Fb != 0.0 && std::copysign(Fa, Fb) == Fa)
    {
      throw std::runtime_error("BrentRootFinder::zeroIn: F(a) and F(b) must have different signs. "
                               "F(a) = " +
                               std::to_string(Fa) + ", F(b) = " + std::to_string(Fb));
    }

    double c = a;
    double Fc = Fa;
    double d = b - a;
    double e = d;

    while (Fb != 0.0)
    {
      // Ensure |F(b)| <= |F(c)|
      if (std::copysign(Fb, Fc) == Fb)
      {
        c = a;
        Fc = Fa;
        d = b - a;
        e = d;
      }

      if (std::fabs(Fc) < std::fabs(Fb))
      {
        a = b;
        b = c;
        c = a;
        Fa = Fb;
        Fb = Fc;
        Fc = Fa;
      }

      // Convergence test
      double xm = 0.5 * (c - b);
      double tol1 = 2.0 * eps * std::fabs(b) + 0.5 * tol;

      if (std::fabs(xm) <= tol1 || Fb == 0.0)
        break;

      // Decide between bisection and interpolation
      if (std::fabs(e) < tol1 || std::fabs(Fa) <= std::fabs(Fb))
      {
        // Bisection
        d = xm;
        e = d;
      }
      else
      {
        double s = Fb / Fa;
        double p, q;

        if (a != c)
        {
          // Inverse quadratic interpolation
          double q1 = Fa / Fc;
          double r = Fb / Fc;
          p = s * (2.0 * xm * q1 * (q1 - r) - (b - a) * (r - 1.0));
          q = (q1 - 1.0) * (r - 1.0) * (s - 1.0);
        }
        else
        {
          // Linear interpolation (secant method)
          p = 2.0 * xm * s;
          q = 1.0 - s;
        }

        // Adjust signs
        if (p > 0.0)
          q = -q;
        else
          p = -p;

        // Accept interpolation?
        if (2.0 * p < 3.0 * xm * q - std::fabs(tol1 * q) && p < std::fabs(0.5 * e * q))
        {
          e = d;
          d = p / q;
        }
        else
        {
          // Interpolation failed, use bisection
          d = xm;
          e = d;
        }
      }

      // Move to new iterate
      a = b;
      Fa = Fb;

      if (std::fabs(d) > tol1)
        b += d;
      else
        b += std::copysign(tol1, xm);

      Fb = F(b);
    }

    return b;
  }

  /**
   * Solve for slip rate V from the quasi-dynamic traction balance equation:
   *   τ_qs = σn * f(V, θ) + η * V
   *
   * where f(V, θ) is the regularized friction coefficient.
   *
   * @param tau_qs Quasi-static shear traction from elasticity (Pa)
   * @param sigma_n Normal stress, positive in compression (Pa)
   * @param theta State variable (s)
   * @param eta Radiation damping coefficient (Pa·s/m)
   * @param a Direct effect parameter
   * @param b Evolution effect parameter (unused here, but needed for f)
   * @param Dc Critical slip distance (m)
   * @param f0 Reference friction coefficient
   * @param V0 Reference slip velocity (m/s)
   * @return Slip rate V (m/s)
   */
  static double solveSlipRate(double tau_qs,
                              double sigma_n,
                              double theta,
                              double eta,
                              double a,
                              double b,
                              double Dc,
                              double f0,
                              double V0)
  {
    // For tension (sigma_n <= 0), slip rate is simply tau/eta
    if (sigma_n <= 0.0)
      return std::fabs(tau_qs) / eta;

    // Define the residual function: R(V) = τ_qs - σn * f(V, θ) - η * V
    auto residual = [&](double V) {
      // Regularized friction coefficient
      double arg = V / (2.0 * V0) * std::exp((f0 + b * std::log(V0 * theta / Dc)) / a);
      double f = a * std::asinh(arg);
      return tau_qs - sigma_n * f - eta * V;
    };

    // Bracket for V: [0, τ_qs/η]
    // At V = 0: R(0) = τ_qs - σn * f(0, θ) ≈ τ_qs (since f(0) → 0 for regularized law)
    // At V = τ_qs/η: R = τ_qs - σn * f(V) - τ_qs = -σn * f(V) < 0
    double V_lower = 1e-20; // Small positive to avoid log(0)
    double V_upper = std::fabs(tau_qs) / eta;

    // Ensure valid bracket
    if (V_upper <= V_lower)
      V_upper = V_lower + 1.0;

    // Check signs at boundaries
    double R_lower = residual(V_lower);
    double R_upper = residual(V_upper);

    // If same sign, expand the bracket
    if (R_lower * R_upper > 0)
    {
      // Try larger upper bound
      V_upper = std::max(V_upper * 10.0, 1.0);
      R_upper = residual(V_upper);

      if (R_lower * R_upper > 0)
      {
        // Fall back to a simple estimate
        // τ_qs ≈ f * σn + η * V, ignoring friction: V ≈ τ_qs / η
        return std::fabs(tau_qs) / (eta + sigma_n * a / V0);
      }
    }

    // Find the root
    return zeroIn(V_lower, V_upper, residual);
  }
};
