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

class LocalizingCoefDiffusion : public Kernel
{
public:
  static InputParameters validParams();

  LocalizingCoefDiffusion(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  Real computeinteractionfunc();

  const Real & _coef;
  const Real & _R;
  const Real & _eta;

  const MaterialProperty<Real> & _d;

};