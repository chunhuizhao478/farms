//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CoupledElkLocalEqstrainForce.h"

registerMooseObject("farmsApp", CoupledElkLocalEqstrainForce);

InputParameters
CoupledElkLocalEqstrainForce::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Kernel for implement local equivalent strain force");
  params.addRequiredCoupledVar(
      "eqstrain_local",
      "The local equivalent strain used in the damage evolution law");
  return params;
}

CoupledElkLocalEqstrainForce::CoupledElkLocalEqstrainForce(const InputParameters & parameters)
  : Kernel(parameters),
    _eqstrain_local(coupledValue("eqstrain_local"))
{
}

Real
CoupledElkLocalEqstrainForce::computeQpResidual()
{
  return -1.0 * _test[_i][_qp] * _eqstrain_local[_qp];
}

Real
CoupledElkLocalEqstrainForce::computeQpJacobian()
{
  return 0.0;
}