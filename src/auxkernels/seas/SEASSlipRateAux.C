//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SEASSlipRateAux.h"
#include "BrentRootFinder.h"
#include <cmath>

registerMooseObject("farmsApp", SEASSlipRateAux);

InputParameters
SEASSlipRateAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription(
      "Computes slip rate V from traction balance using Brent's method. "
      "Solves: τ = σn * f(V, θ) + η * V for V.");

  params.addRequiredCoupledVar("traction", "The elastic traction variable");
  params.addRequiredCoupledVar("state_variable", "The state variable θ");

  // Rate-state parameters
  params.addRequiredParam<Real>("a", "Direct effect parameter a");
  params.addRequiredParam<Real>("b", "Evolution effect parameter b");
  params.addRequiredParam<Real>("Dc", "Critical slip distance Dc (m)");
  params.addParam<Real>("f0", 0.6, "Reference friction coefficient");
  params.addParam<Real>("V0", 1e-6, "Reference slip velocity (m/s)");
  params.addRequiredParam<Real>("sigma_n", "Normal stress, positive in compression (Pa)");

  // Radiation damping
  params.addRequiredParam<Real>("eta", "Radiation damping coefficient η = μ/(2*cs) (Pa·s/m)");

  // Regularization
  params.addParam<Real>("V_min", 1e-20, "Minimum slip rate for regularization (m/s)");

  return params;
}

SEASSlipRateAux::SEASSlipRateAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _traction(coupledValue("traction")),
    _state_variable(coupledValue("state_variable")),
    _a(getParam<Real>("a")),
    _b(getParam<Real>("b")),
    _Dc(getParam<Real>("Dc")),
    _f0(getParam<Real>("f0")),
    _V0(getParam<Real>("V0")),
    _sigma_n(getParam<Real>("sigma_n")),
    _eta(getParam<Real>("eta")),
    _V_min(getParam<Real>("V_min"))
{
}

Real
SEASSlipRateAux::computeValue()
{
  Real tau = _traction[_qp];
  Real theta = _state_variable[_qp];

  // Ensure state variable is positive
  if (theta <= 0.0)
    theta = _Dc / _V_min;

  // Solve for slip rate
  return solveSlipRate(tau, theta);
}

Real
SEASSlipRateAux::solveSlipRate(Real tau, Real theta) const
{
  // For tension (sigma_n <= 0), slip rate is simply tau/eta
  if (_sigma_n <= 0.0)
    return std::max(std::fabs(tau) / _eta, _V_min);

  // Handle negative traction (should not happen in mode III, but handle gracefully)
  Real tau_abs = std::fabs(tau);

  // Define the residual function: R(V) = τ - σn * f(V, θ) - η * V
  auto residual = [&](double V) {
    Real f = frictionCoefficient(V, theta);
    return tau_abs - _sigma_n * f - _eta * V;
  };

  // Bracket for V: [V_min, V_max]
  Real V_min = _V_min;
  Real V_max = std::max(tau_abs / _eta, 1.0);

  // Check signs at boundaries to ensure we have a bracket
  Real R_min = residual(V_min);
  Real R_max = residual(V_max);

  // If same sign, expand the bracket
  int max_expand = 10;
  while (R_min * R_max > 0 && max_expand > 0)
  {
    V_max *= 10.0;
    R_max = residual(V_max);
    max_expand--;
  }

  // If still no valid bracket, use linear approximation
  if (R_min * R_max > 0)
  {
    // τ ≈ (a * σn / V0) * V + η * V for small V
    return std::max(tau_abs / (_a * _sigma_n / _V0 + _eta), _V_min);
  }

  // Ensure proper bracket ordering (R should go from + to -)
  if (R_min < R_max)
  {
    std::swap(V_min, V_max);
  }

  // Find the root using Brent's method
  try
  {
    return BrentRootFinder::zeroIn(V_min, V_max, residual);
  }
  catch (const std::exception &)
  {
    // Fallback to linear approximation
    return std::max(tau_abs / (_a * _sigma_n / _V0 + _eta), _V_min);
  }
}

Real
SEASSlipRateAux::frictionCoefficient(Real V, Real theta) const
{
  // Regularized friction coefficient (Tandem paper Eq. 7):
  // f(V, θ) = a * asinh[V/(2*V0) * exp((f0 + b*ln(V0*θ/Dc))/a)]

  // Ensure positive arguments for log
  V = std::max(V, _V_min);
  theta = std::max(theta, _Dc / _V0);

  Real arg = V / (2.0 * _V0) * std::exp((_f0 + _b * std::log(_V0 * theta / _Dc)) / _a);
  return _a * std::asinh(arg);
}
