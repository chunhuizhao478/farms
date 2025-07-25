//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Kernel.h"

/**
 * This kernel implements a simple reaction term with a coupled local equivalent strain
 */
class CoupledReaction : public Kernel
{
public:
  static InputParameters validParams();
  CoupledReaction(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  /// Reaction rate
  const Real & _rate;
  
  /// Coupled local equivalent strain
  const VariableValue & _eqstrain_local;
  /// Length scale for gradient activity parameter
  const Real _length_scale;
  /// Equivalent strain at which the gradient activity starts
  const Real _kappa_i;  
  /// Minimum value of the gradient activity parameter for the equivalent strain
  const Real _c0;
};