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
 * SEASSlipAux integrates the slip rate to update the accumulated slip.
 *
 * Simple forward Euler:
 *   S_new = S_old + dt * V
 *
 * This is part of the staggered SEAS solver where friction is updated explicitly.
 */
class SEASSlipAux : public AuxKernel
{
public:
  static InputParameters validParams();

  SEASSlipAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  /// Slip rate (from SEASSlipRateAux)
  const VariableValue & _slip_rate;

  /// Old slip value
  const VariableValue & _slip_old;
};
