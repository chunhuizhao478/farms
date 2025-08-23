//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LocalizingCoefDiffusion.h"

registerMooseObject("farmsApp", LocalizingCoefDiffusion);

InputParameters
LocalizingCoefDiffusion::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addCustomTypeParam("coef", 0.0, "CoefficientType", "The coefficient of diffusion");
  params.addPrivateParam<Real>("_test_private_param", 12345);
  params.addParam<Real>("non_controllable", "A parameter we cannot control.");

  params.declareControllable("coef");

  params.addParam<Real>("R", 0.005, "The internal length scale");
  params.addParam<Real>("eta", 5, "The internal length scale");

  return params;
}

LocalizingCoefDiffusion::LocalizingCoefDiffusion(const InputParameters & parameters)
  : Kernel(parameters), 
  _coef(getParam<Real>("coef")),
  _R(getParam<Real>("R")),
  _eta(getParam<Real>("eta")),
  _d(getMaterialProperty<Real>("crack_damage"))
{
}

Real
LocalizingCoefDiffusion::computeQpResidual()
{
  Real g = computeinteractionfunc();
  return _coef * _grad_test[_i][_qp] * _grad_u[_qp] * g;
}

Real
LocalizingCoefDiffusion::computeQpJacobian()
{
  Real g = computeinteractionfunc();
  return _coef * _grad_test[_i][_qp] * _grad_phi[_j][_qp] * g;
}

Real
LocalizingCoefDiffusion::computeinteractionfunc()
{
   return ((1 - _R) * exp(-_eta * _d[_qp]) + _R - exp(-_eta)) / (1 - exp(-_eta));
}