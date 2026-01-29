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
 * SEASSlipRateVarAAux computes slip rate V from traction balance using Brent's method.
 * This version supports spatially-varying rate-state parameter 'a' via a coupled variable.
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
 * Key difference from SEASSlipRateAux: 'a' is read from a coupled AuxVariable
 * allowing spatially-varying a(z) as required by BP2 benchmark.
 *
 * Backslip Loading (Tandem-style):
 *   For depths z > backslip_depth, prescribes V = Vp (plate rate)
 *   instead of solving the friction balance. This creates the correct
 *   stress loading pattern for half-space SEAS simulations.
 */
class SEASSlipRateVarAAux : public AuxKernel
{
public:
  static InputParameters validParams();

  SEASSlipRateVarAAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  /// Solve for slip rate from traction balance
  Real solveSlipRate(Real tau, Real theta, Real a_local) const;

  /// Compute regularized friction coefficient
  Real frictionCoefficient(Real V, Real theta, Real a_local) const;

  /// Elastic traction (from SEASTractionAux)
  const VariableValue & _traction;

  /// State variable (from SEASStateAux)
  const VariableValue & _state_variable;

  /// Rate-state parameter 'a' (spatially-varying, from AuxVariable)
  const VariableValue & _a_var;

  // Rate-state parameters (constant)
  const Real _b;        ///< Evolution effect parameter
  const Real _Dc;       ///< Critical slip distance (m)
  const Real _f0;       ///< Reference friction coefficient
  const Real _V0;       ///< Reference slip velocity (m/s)
  const Real _sigma_n;  ///< Normal stress (Pa), positive in compression

  /// Radiation damping coefficient η = μ / (2 * cs)
  const Real _eta;

  /// Minimum slip rate (regularization)
  const Real _V_min;

  /// Maximum slip rate (stability limit)
  const Real _V_max;

  /// Pre-stress (initial shear stress)
  const Real _tau_pre;

  // Backslip loading parameters (Tandem-style)
  /// Enable backslip loading
  const bool _use_backslip;

  /// Plate rate for backslip (m/s)
  const Real _Vp;

  /// Depth below which backslip is prescribed (m)
  const Real _backslip_depth;

  /// Coordinate direction for depth (0=x, 1=y, 2=z)
  const unsigned int _depth_dir;
};
