//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PorePressureTimeDerivative.h"

registerMooseObject("farmsApp", PorePressureTimeDerivative);

InputParameters
PorePressureTimeDerivative::validParams()
{
  InputParameters params = TimeKernel::validParams();
  params.addClassDescription("Time derivative kernel for pore pressure in poromechanical problems");
  
  return params;
}

PorePressureTimeDerivative::PorePressureTimeDerivative(const InputParameters & parameters)
  : TimeKernel(parameters),
    _u_old(valueOld()),
    _Biot_modulus_eff(getMaterialProperty<Real>("Biot_modulus_effective"))
{
}

Real
PorePressureTimeDerivative::computeQpResidual()
{
  return _test[_i][_qp] * (1.0 / _Biot_modulus_eff[_qp]) *  (_u[_qp] - _u_old[_qp]) / _dt ;
}

Real
PorePressureTimeDerivative::computeQpJacobian()
{
  return _test[_i][_qp] * (1.0 / _Biot_modulus_eff[_qp])  * _phi[_j][_qp];
}