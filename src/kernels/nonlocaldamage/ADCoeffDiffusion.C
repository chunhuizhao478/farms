//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADCoefDiffusion.h"

registerMooseObject("farmsApp", ADCoefDiffusion);

InputParameters
ADCoefDiffusion::validParams()
{
  InputParameters params = ADKernel::validParams();
  params.addParam<Real>("coef", 0.0, "Diffusion coefficient");
  params.addParam<FunctionName>("function",
                                "If provided, the diffusion coefficient will be coef + "
                                "this function.  This is useful for temporally or "
                                "spatially varying diffusivities");
  params.addClassDescription("Kernel for diffusion with diffusivity = coef + function");
  return params;
}

ADCoefDiffusion::ADCoefDiffusion(const InputParameters & parameters)
  : ADKernel(parameters),
    _coef(getParam<Real>("coef")),
    _func(parameters.isParamValid("function") ? &getFunction("function") : NULL)
{
}

ADReal
ADCoefDiffusion::computeQpResidual()
{
  ADReal diffusivity = _coef;

  if (_func)
    diffusivity += _func->value(_t, _q_point[_qp]);

  return diffusivity * _grad_test[_i][_qp] * _grad_u[_qp];
}