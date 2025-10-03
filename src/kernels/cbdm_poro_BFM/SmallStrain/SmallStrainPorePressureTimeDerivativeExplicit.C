//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SmallStrainPorePressureTimeDerivativeExplicit.h"

registerMooseObject("farmsApp", SmallStrainPorePressureTimeDerivativeExplicit);

InputParameters
SmallStrainPorePressureTimeDerivativeExplicit::validParams()
{
  InputParameters params = TimeKernel::validParams();
  params.addClassDescription("Time derivative kernel for pore pressure in poromechanical problems");
  
  return params;
}

SmallStrainPorePressureTimeDerivativeExplicit::SmallStrainPorePressureTimeDerivativeExplicit(const InputParameters & parameters)
  : TimeKernel(parameters),
    _u_older(valueOlder()),
    _Biot_modulus_eff(getMaterialProperty<Real>("Biot_modulus_effective"))
{
}

Real
SmallStrainPorePressureTimeDerivativeExplicit::computeQpResidual()
{
  return _test[_i][_qp] * (1.0 / _Biot_modulus_eff[_qp]) *  (_u[_qp] - _u_older[_qp]) / _dt / 2;
}

Real
SmallStrainPorePressureTimeDerivativeExplicit::computeQpJacobian()
{
  return _test[_i][_qp] * (1.0 / _Biot_modulus_eff[_qp])  * _phi[_j][_qp] / _dt / 2;
}