//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FarmsLocalXiForce.h"

registerMooseObject("farmsApp", FarmsLocalXiForce);

InputParameters
FarmsLocalXiForce::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Kernel for implement local equivalent strain force");
  return params;
}

FarmsLocalXiForce::FarmsLocalXiForce(const InputParameters & parameters)
  : Kernel(parameters),
    _eqstrain_local(getMaterialProperty<Real>("strain_invariant_ratio"))
{
}

Real
FarmsLocalXiForce::computeQpResidual()
{
  return -1.0 * _test[_i][_qp] * _eqstrain_local[_qp];
}

Real
FarmsLocalXiForce::computeQpJacobian()
{
  return 0.0;
}