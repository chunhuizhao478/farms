//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SEASStateAux.h"
#include <cmath>

registerMooseObject("farmsApp", SEASStateAux);

InputParameters
SEASStateAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription(
      "Updates state variable θ using the aging law with backward Euler. "
      "dθ/dt = 1 - V*θ/Dc => θ_new = (θ_old + dt) / (1 + V*dt/Dc)");

  params.addRequiredCoupledVar("slip_rate", "The slip rate variable V");
  params.addRequiredParam<Real>("Dc", "Critical slip distance Dc (m)");

  return params;
}

SEASStateAux::SEASStateAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _slip_rate(coupledValue("slip_rate")),
    _state_old(_var.slnOld()),
    _Dc(getParam<Real>("Dc"))
{
}

Real
SEASStateAux::computeValue()
{
  Real V = std::fabs(_slip_rate[_qp]);
  Real theta_old = _state_old[_qp];

  // Ensure positive state
  if (theta_old <= 0.0)
    theta_old = _Dc / std::max(V, 1e-20);

  // Backward Euler for aging law: dθ/dt = 1 - V*θ/Dc
  // θ_new = (θ_old + dt) / (1 + V*dt/Dc)
  if (_dt > 0)
  {
    return (theta_old + _dt) / (1.0 + V * _dt / _Dc);
  }
  else
  {
    // Initial condition
    return theta_old;
  }
}
