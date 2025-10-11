//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FarmsVelocityMagnitudeAux.h"

#include <cmath>

registerMooseObject("farmsApp", FarmsVelocityMagnitudeAux);

InputParameters
FarmsVelocityMagnitudeAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("Computes velocity magnitude from component auxiliary variables.");
  params.addRequiredCoupledVar("vel_x", "Velocity component in x direction.");
  params.addRequiredCoupledVar("vel_y", "Velocity component in y direction.");
  params.addCoupledVar("vel_z", "Velocity component in z direction (optional).");
  return params;
}

FarmsVelocityMagnitudeAux::FarmsVelocityMagnitudeAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _vel_x(coupledValue("vel_x")),
    _vel_y(coupledValue("vel_y")),
    _vel_z(isCoupled("vel_z") ? &coupledValue("vel_z") : nullptr),
    _has_z(isCoupled("vel_z"))
{
}

Real
FarmsVelocityMagnitudeAux::computeValue()
{
  const Real vx = _vel_x[_qp];
  const Real vy = _vel_y[_qp];
  const Real vz = _has_z ? (*_vel_z)[_qp] : 0.0;
  return std::sqrt(vx * vx + vy * vy + vz * vz);
}
