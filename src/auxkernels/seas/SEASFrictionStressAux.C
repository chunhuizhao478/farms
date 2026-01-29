//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SEASFrictionStressAux.h"
#include <cmath>

registerMooseObject("farmsApp", SEASFrictionStressAux);

InputParameters
SEASFrictionStressAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription(
      "Computes shear stress from friction law: τ = σn * f(V, θ) + η * V. "
      "This gives the physically consistent stress for SEAS simulations.");

  params.addRequiredCoupledVar("slip_rate", "The slip rate variable V");
  params.addRequiredCoupledVar("state_variable", "The state variable θ");
  params.addRequiredCoupledVar("a_var", "Spatially-varying rate-state parameter 'a'");

  // Rate-state parameters
  params.addRequiredParam<Real>("b", "Evolution effect parameter b");
  params.addParam<Real>("f0", 0.6, "Reference friction coefficient");
  params.addParam<Real>("V0", 1e-6, "Reference slip velocity (m/s)");
  params.addParam<Real>("Dc", 0.004, "Critical slip distance Dc (m)");
  params.addRequiredParam<Real>("sigma_n", "Normal stress, positive in compression (Pa)");
  params.addRequiredParam<Real>("eta", "Radiation damping coefficient η = μ/(2*cs) (Pa·s/m)");

  return params;
}

SEASFrictionStressAux::SEASFrictionStressAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _slip_rate(coupledValue("slip_rate")),
    _state_variable(coupledValue("state_variable")),
    _a_var(coupledValue("a_var")),
    _b(getParam<Real>("b")),
    _f0(getParam<Real>("f0")),
    _V0(getParam<Real>("V0")),
    _Dc(getParam<Real>("Dc")),
    _sigma_n(getParam<Real>("sigma_n")),
    _eta(getParam<Real>("eta"))
{
}

Real
SEASFrictionStressAux::computeValue()
{
  Real V = _slip_rate[_qp];
  Real theta = _state_variable[_qp];
  Real a = _a_var[_qp];

  // Ensure positive values
  V = std::max(V, 1e-20);
  theta = std::max(theta, _Dc / _V0);
  if (a <= 0.0)
    a = 0.025;

  // Compute friction coefficient using regularized form (Tandem Eq. 7):
  // f(V, θ) = a * asinh[V/(2*V0) * exp((f0 + b*ln(V0*θ/Dc))/a)]
  Real arg = V / (2.0 * _V0) * std::exp((_f0 + _b * std::log(_V0 * theta / _Dc)) / a);
  Real f = a * std::asinh(arg);

  // Shear stress from friction law: τ = σn * f(V, θ) + η * V
  return _sigma_n * f + _eta * V;
}
