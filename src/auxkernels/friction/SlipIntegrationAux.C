//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SlipIntegrationAux.h"

registerMooseObject("farmsApp", SlipIntegrationAux);

InputParameters
SlipIntegrationAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription(
      "Integrates slip rate to update accumulated slip: S_new = S_old + dt * V. "
      "Designed for the Friction SubApp in MultiApp SEAS architecture. "
      "Uses coupledValueOld for proper state preservation in MultiApp context.");

  params.addRequiredCoupledVar("slip_rate", "The slip rate variable V (m/s)");
  params.addRequiredCoupledVar("slip_old_var", "The slip variable to get old value from (couple to self)");

  return params;
}

SlipIntegrationAux::SlipIntegrationAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _slip_rate(coupledValue("slip_rate")),
    _slip_old(coupledValueOld("slip_old_var"))
{
}

Real
SlipIntegrationAux::computeValue()
{
  Real V = _slip_rate[_qp];
  Real slip_old = _slip_old[_qp];

  // Forward Euler: S_new = S_old + dt * V
  if (_dt > 0)
  {
    return slip_old + _dt * V;
  }
  else
  {
    // Initial condition or zero dt
    return slip_old;
  }
}
