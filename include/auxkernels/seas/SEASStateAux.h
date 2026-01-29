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
 * SEASStateAux updates the state variable θ using the aging law.
 *
 * Aging law:
 *   dθ/dt = 1 - V*θ/Dc
 *
 * Discretized with backward Euler:
 *   θ_new = (θ_old + dt) / (1 + V*dt/Dc)
 *
 * This provides unconditional stability for the state evolution.
 *
 * This is part of the staggered SEAS solver where friction is updated explicitly.
 */
class SEASStateAux : public AuxKernel
{
public:
  static InputParameters validParams();

  SEASStateAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  /// Slip rate (from SEASSlipRateAux)
  const VariableValue & _slip_rate;

  /// Old state variable value
  const VariableValue & _state_old;

  /// Critical slip distance
  const Real _Dc;
};
