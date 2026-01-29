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
 * SEASFrictionStressAux computes the shear stress on the fault from the
 * friction law: τ = σn * f(V, θ) + η * V
 *
 * This gives the physically consistent stress that satisfies the traction
 * balance equation in SEAS simulations. The friction coefficient f(V, θ)
 * uses the regularized rate-and-state formulation.
 */
class SEASFrictionStressAux : public AuxKernel
{
public:
  static InputParameters validParams();

  SEASFrictionStressAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  /// Slip rate variable
  const VariableValue & _slip_rate;

  /// State variable
  const VariableValue & _state_variable;

  /// Spatially-varying a parameter
  const VariableValue & _a_var;

  /// Rate-state parameter b
  const Real _b;

  /// Reference friction coefficient
  const Real _f0;

  /// Reference slip velocity
  const Real _V0;

  /// Critical slip distance
  const Real _Dc;

  /// Normal stress (positive in compression)
  const Real _sigma_n;

  /// Radiation damping coefficient
  const Real _eta;
};
