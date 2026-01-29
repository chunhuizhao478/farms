//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "InterfaceMaterial.h"

/**
 * DGRateStateFrictionMaterial implements rate-and-state friction
 * for DG SEAS simulations on fault interfaces.
 *
 * This material uses the Tandem-style decoupled approach:
 * 1. Compute elastic traction from displacement gradients
 * 2. Solve for slip rate V from traction balance using Brent's method
 * 3. Update state variable using aging law
 *
 * Regularized friction coefficient (Tandem paper Eq. 7):
 *   f(V,θ) = a * asinh[V/(2*V0) * exp((f0 + b*ln(V0*θ/Dc))/a)]
 *
 * Aging law for state evolution:
 *   dθ/dt = 1 - V*θ/Dc
 *
 * Traction balance (quasi-dynamic):
 *   τ_elastic = σn * f(V,θ) + η * V
 *
 * This material computes:
 * - fault_traction: shear traction from friction law (for interface kernel)
 * - slip_rate: slip velocity V (solved from traction balance)
 * - state_variable: state variable θ
 * - dtraction_dslip: derivative for Jacobian
 */
class DGRateStateFrictionMaterial : public InterfaceMaterial
{
public:
  static InputParameters validParams();

  DGRateStateFrictionMaterial(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties() override;
  virtual void computeQpProperties() override;

  /// Solve for slip rate V given elastic traction and state
  Real solveSlipRate(Real tau_elastic, Real theta) const;

  /// Compute friction coefficient
  Real frictionCoefficient(Real V, Real theta) const;

  /// Compute derivative of friction coefficient w.r.t. slip rate
  Real dFrictionDSlipRate(Real V, Real theta) const;

  /// Compute elastic traction from displacement gradients
  Real computeElasticTraction() const;

  // Input variables
  /// Displacement on element side
  const VariableValue & _u;
  /// Displacement on neighbor side
  const VariableValue & _u_neighbor;
  /// Displacement gradient on element side
  const VariableGradient & _grad_u;
  /// Displacement gradient on neighbor side
  const VariableGradient & _grad_u_neighbor;

  // Normal vector
  const MooseArray<Point> & _normals;

  // Rate-state friction parameters
  const Real _a;        ///< Direct effect parameter
  const Real _b;        ///< Evolution effect parameter
  const Real _Dc;       ///< Critical slip distance (m)
  const Real _f0;       ///< Reference friction coefficient
  const Real _V0;       ///< Reference slip velocity (m/s)
  const Real _sigma_n;  ///< Normal stress (Pa), positive in compression

  // Material properties
  const Real _shear_modulus;
  const Real _density;

  // Initial slip rate for initialization
  const Real _initial_slip_rate;

  // Output material properties
  MaterialProperty<Real> & _fault_traction;
  MaterialProperty<Real> & _slip;
  MaterialProperty<Real> & _slip_rate;
  MaterialProperty<Real> & _state_variable;
  MaterialProperty<Real> & _dtraction_dslip;
  MaterialProperty<Real> & _friction_coefficient;

  // Old values for time integration
  const MaterialProperty<Real> & _state_variable_old;
  const MaterialProperty<Real> & _slip_old;
  const MaterialProperty<Real> & _slip_rate_old;

  // Radiation damping coefficient
  Real _eta;

  // Pre-stress (optional, for initialization)
  const Real _tau_pre;
};
