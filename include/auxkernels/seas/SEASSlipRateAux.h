//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "AuxKernel.h"

/**
 * SEASSlipRateAux computes the slip rate V from the traction balance equation
 * using Brent's method for robust root finding.
 *
 * Solves: τ = σn * f(V, θ) + η * V
 *
 * where:
 *   τ = elastic traction (from displacement field)
 *   σn = normal stress (positive in compression)
 *   f(V, θ) = regularized friction coefficient
 *   η = radiation damping coefficient = μ / (2 * cs)
 *   V = slip rate (unknown)
 *   θ = state variable
 *
 * Regularized friction (Tandem Eq. 7):
 *   f(V, θ) = a * asinh[V/(2*V0) * exp((f0 + b*ln(V0*θ/Dc))/a)]
 *
 * This is part of the staggered SEAS solver where friction is updated explicitly.
 */
class SEASSlipRateAux : public AuxKernel
{
public:
  static InputParameters validParams();

  SEASSlipRateAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  /// Solve for slip rate from traction balance
  Real solveSlipRate(Real tau, Real theta) const;

  /// Compute regularized friction coefficient
  Real frictionCoefficient(Real V, Real theta) const;

  /// Elastic traction (from SEASTractionAux)
  const VariableValue & _traction;

  /// State variable (from SEASStateAux)
  const VariableValue & _state_variable;

  // Rate-state parameters
  const Real _a;        ///< Direct effect parameter
  const Real _b;        ///< Evolution effect parameter
  const Real _Dc;       ///< Critical slip distance (m)
  const Real _f0;       ///< Reference friction coefficient
  const Real _V0;       ///< Reference slip velocity (m/s)
  const Real _sigma_n;  ///< Normal stress (Pa), positive in compression

  /// Radiation damping coefficient η = μ / (2 * cs)
  const Real _eta;

  /// Minimum slip rate (regularization)
  const Real _V_min;
};
