//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ADKernel.h"
#include "Function.h"

class ADCoefDiffusion : public ADKernel
{
public:
  static InputParameters validParams();

  ADCoefDiffusion(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual();

private:
  const ADReal _coef;
  const Function * const _func;
};