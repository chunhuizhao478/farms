//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FarmsVelocityAux.h"

registerMooseObject("farmsApp", FarmsVelocityAux);

InputParameters
FarmsVelocityAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("Populates a velocity auxiliary variable with the time derivative of a displacement.");
  params.addRequiredCoupledVar("displacement", "Displacement variable whose time derivative defines the velocity.");
  return params;
}

FarmsVelocityAux::FarmsVelocityAux(const InputParameters & parameters)
  : AuxKernel(parameters), _disp_dot(coupledDot("displacement"))
{
}

Real
FarmsVelocityAux::computeValue()
{
  return _disp_dot[_qp];
}
