//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SEASAdaptiveDT.h"
#include <cmath>
#include <algorithm>

registerMooseObject("farmsApp", SEASAdaptiveDT);

InputParameters
SEASAdaptiveDT::validParams()
{
  InputParameters params = TimeStepper::validParams();
  params.addClassDescription(
      "Adaptive time stepping for SEAS simulations. Adjusts dt based on "
      "maximum slip rate: small dt during seismic, large dt during aseismic periods.");

  params.addRequiredParam<PostprocessorName>(
      "max_slip_rate_pp", "Postprocessor providing maximum slip rate on the fault");

  params.addRequiredParam<Real>("Dc", "Critical slip distance (m)");
  params.addParam<Real>("C", 0.5, "Safety factor for time step calculation (dt = C * Dc / V)");

  params.addParam<Real>("dt_min", 1e-3, "Minimum allowed time step (s)");
  params.addParam<Real>("dt_max", 3.15576e7, "Maximum allowed time step (s), default ~1 year");

  params.addParam<Real>("V_seismic", 1e-3,
                        "Threshold slip rate (m/s) for seismic vs aseismic classification");
  params.addParam<Real>("dt_seismic", 0.1,
                        "Target time step during seismic phase (s)");

  params.addParam<Real>("initial_dt", 1.0, "Initial time step (s)");
  params.addParam<Real>("growth_factor", 1.2,
                        "Maximum growth factor for dt between time steps during aseismic phase");

  return params;
}

SEASAdaptiveDT::SEASAdaptiveDT(const InputParameters & parameters)
  : TimeStepper(parameters),
    PostprocessorInterface(this),
    _max_slip_rate(getPostprocessorValue("max_slip_rate_pp")),
    _Dc(getParam<Real>("Dc")),
    _C(getParam<Real>("C")),
    _dt_min(getParam<Real>("dt_min")),
    _dt_max(getParam<Real>("dt_max")),
    _V_seismic(getParam<Real>("V_seismic")),
    _dt_seismic(getParam<Real>("dt_seismic")),
    _initial_dt(getParam<Real>("initial_dt")),
    _growth_factor(getParam<Real>("growth_factor"))
{
}

Real
SEASAdaptiveDT::computeInitialDT()
{
  return _initial_dt;
}

Real
SEASAdaptiveDT::computeDT()
{
  // Get current maximum slip rate (ensure positive value)
  Real V_max = std::abs(_max_slip_rate);

  // Regularize to avoid division by zero
  if (V_max < 1e-20)
    V_max = 1e-20;

  Real dt_new;
  Real dt_old = getCurrentDT();
  Real dt_slip = _C * _Dc / V_max;

  // Check if we're in seismic or aseismic phase
  if (V_max >= _V_seismic)
  {
    // Seismic phase: use small, controlled time step
    // dt = C * Dc / V_max, but capped at dt_seismic
    dt_new = std::min(dt_slip, _dt_seismic);
  }
  else
  {
    // Aseismic phase: compute dt based on slip rate
    dt_new = dt_slip;

    // Limit growth factor to ensure smooth time stepping
    if (dt_old > 0)
    {
      dt_new = std::min(dt_new, _growth_factor * dt_old);
    }
  }

  // Apply bounds
  dt_new = std::max(dt_new, _dt_min);
  dt_new = std::min(dt_new, _dt_max);

  // Debug output
  _console << "SEASAdaptiveDT: V_max=" << V_max << " dt_slip=" << dt_slip
           << " dt_old=" << dt_old << " growth_limit=" << _growth_factor * dt_old
           << " dt_new=" << dt_new << std::endl;

  return dt_new;
}
