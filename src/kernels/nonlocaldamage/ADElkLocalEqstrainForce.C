//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADElkLocalEqstrainForce.h"

registerMooseObject("farmsApp", ADElkLocalEqstrainForce);

InputParameters
ADElkLocalEqstrainForce::validParams()
{
  InputParameters params = ADKernel::validParams();
  params.addClassDescription("Kernel for implement local equivalent strain force");
  return params;
}

ADElkLocalEqstrainForce::ADElkLocalEqstrainForce(const InputParameters & parameters)
  : ADKernel(parameters),
    _eqstrain_local(getADMaterialProperty<Real>("eqstrain_local"))
{
}

ADReal
ADElkLocalEqstrainForce::computeQpResidual()
{
  return -1.0 * _test[_i][_qp] * _eqstrain_local[_qp];
}