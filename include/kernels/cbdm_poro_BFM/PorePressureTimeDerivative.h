//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "TimeKernel.h"

class PorePressureTimeDerivative : public TimeKernel
{
public:
  static InputParameters validParams();

  PorePressureTimeDerivative(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  /// old value of porepressure
  const VariableValue & _u_old;

  /// Biot modulus
  const MaterialProperty<Real> & _Biot_modulus_eff;
};