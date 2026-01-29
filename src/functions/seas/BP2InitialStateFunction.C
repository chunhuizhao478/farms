//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "BP2InitialStateFunction.h"
#include <cmath>

registerMooseObject("farmsApp", BP2InitialStateFunction);

InputParameters
BP2InitialStateFunction::validParams()
{
  InputParameters params = Function::validParams();
  params.addClassDescription(
      "Computes initial state variable θ(z,0) for BP2 benchmark using Eq. 12.");

  params.addRequiredParam<Real>("Dc", "Critical slip distance (m)");
  params.addRequiredParam<Real>("V0", "Reference slip rate (m/s)");
  params.addRequiredParam<Real>("Vinit", "Initial slip rate (m/s)");
  params.addRequiredParam<Real>("b", "Rate-state parameter b");
  params.addRequiredParam<Real>("f0", "Reference friction coefficient");
  params.addRequiredParam<Real>("tau0", "Pre-stress / initial shear stress (Pa)");
  params.addRequiredParam<Real>("sigma_n", "Normal stress (Pa)");
  params.addRequiredParam<Real>("eta", "Radiation damping coefficient (Pa·s/m)");
  params.addRequiredParam<Real>("a0", "Rate-state a in VW region");
  params.addRequiredParam<Real>("amax", "Rate-state a in VS region");
  params.addRequiredParam<Real>("H", "Depth of VW region (m)");
  params.addRequiredParam<Real>("h", "Width of VW-VS transition zone (m)");
  params.addParam<unsigned int>("depth_direction", 1, "Coordinate direction for depth (0=x, 1=y, 2=z)");

  return params;
}

BP2InitialStateFunction::BP2InitialStateFunction(const InputParameters & parameters)
  : Function(parameters),
    _Dc(getParam<Real>("Dc")),
    _V0(getParam<Real>("V0")),
    _Vinit(getParam<Real>("Vinit")),
    _b(getParam<Real>("b")),
    _f0(getParam<Real>("f0")),
    _tau0(getParam<Real>("tau0")),
    _sigma_n(getParam<Real>("sigma_n")),
    _eta(getParam<Real>("eta")),
    _a0(getParam<Real>("a0")),
    _amax(getParam<Real>("amax")),
    _H(getParam<Real>("H")),
    _h(getParam<Real>("h")),
    _depth_dir(getParam<unsigned int>("depth_direction"))
{
}

Real
BP2InitialStateFunction::computeA(Real z) const
{
  if (z < _H)
    return _a0;
  else if (z < _H + _h)
    return _a0 + (_amax - _a0) * (z - _H) / _h;
  else
    return _amax;
}

Real
BP2InitialStateFunction::value(Real /*t*/, const Point & p) const
{
  // Get depth from the appropriate coordinate
  // Note: depth coordinate may be negative (y < 0), so use abs()
  Real z = std::fabs(p(_depth_dir));

  // Compute depth-dependent a
  Real a = computeA(z);

  // Compute sinh argument: (τ⁰ - η*Vinit) / (a * σn)
  Real sinh_arg = (_tau0 - _eta * _Vinit) / (a * _sigma_n);

  // Compute the inner term: 2*V0/Vinit * sinh(...)
  Real inner = 2.0 * _V0 / _Vinit * std::sinh(sinh_arg);

  // Compute the exponent: (a/b) * ln(inner) - f0/b
  Real exponent = (a / _b) * std::log(inner) - _f0 / _b;

  // Compute θ = (Dc/V0) * exp(exponent)
  Real theta = (_Dc / _V0) * std::exp(exponent);

  // Debug - no output
  (void)z; (void)a; (void)sinh_arg; (void)inner; (void)exponent; (void)theta;

  return theta;
}
