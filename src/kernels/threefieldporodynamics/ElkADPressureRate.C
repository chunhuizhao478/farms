//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

/**
 * Created by Chunhui Zhao, Jul 3rd, 2024
 * Implement dot(p) / M term
 * dot(p): pressure time derivative (primary variable)
 * M : biot modulus 
 */
#include "ElkADPressureRate.h"

registerMooseObject("farmsApp", ElkADPressureRate);

InputParameters
ElkADPressureRate::validParams()
{
  InputParameters params = ADTimeKernelValue::validParams();
  params.addClassDescription("The time derivative operator with the weak form of $(\\psi_i, "
                             "\\frac{\\partial u_h}{\\partial t})$.");
  params.addParam<bool>("divide_by_biot_modulus", false, "Whether to divide by the Biot modulus");
  return params;
}

ElkADPressureRate::ElkADPressureRate(const InputParameters & parameters)
  : ADTimeKernelValue(parameters),
  _divide_by_biot_modulus(getParam<bool>("divide_by_biot_modulus")),
  _biot_modulus(getADMaterialProperty<Real>("biot_modulus"))
{
}

ADReal
ElkADPressureRate::precomputeQpResidual()
{
  if (_divide_by_biot_modulus)
    return _u_dot[_qp] / _biot_modulus[_qp];
  return _u_dot[_qp];
}