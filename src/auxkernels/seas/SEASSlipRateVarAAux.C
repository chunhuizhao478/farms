//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SEASSlipRateVarAAux.h"
#include "BrentRootFinder.h"
#include <cmath>

registerMooseObject("farmsApp", SEASSlipRateVarAAux);

InputParameters
SEASSlipRateVarAAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription(
      "Computes slip rate V from traction balance using Brent's method. "
      "Supports spatially-varying 'a' parameter via coupled variable. "
      "Solves: τ = σn * f(V, θ) + η * V for V.");

  params.addRequiredCoupledVar("traction", "The elastic traction variable");
  params.addRequiredCoupledVar("state_variable", "The state variable θ");
  params.addRequiredCoupledVar("a_var", "Spatially-varying rate-state parameter 'a'");

  // Rate-state parameters (constant)
  params.addRequiredParam<Real>("b", "Evolution effect parameter b");
  params.addRequiredParam<Real>("Dc", "Critical slip distance Dc (m)");
  params.addParam<Real>("f0", 0.6, "Reference friction coefficient");
  params.addParam<Real>("V0", 1e-6, "Reference slip velocity (m/s)");
  params.addRequiredParam<Real>("sigma_n", "Normal stress, positive in compression (Pa)");

  // Radiation damping
  params.addRequiredParam<Real>("eta", "Radiation damping coefficient η = μ/(2*cs) (Pa·s/m)");

  // Regularization and stability limits
  params.addParam<Real>("V_min", 1e-20, "Minimum slip rate for regularization (m/s)");
  params.addParam<Real>("V_max", 10.0, "Maximum slip rate for stability (m/s)");

  // Pre-stress
  params.addParam<Real>("tau_pre", 0.0, "Initial/pre-stress shear traction (Pa)");

  // Backslip loading parameters (Tandem-style)
  params.addParam<bool>("use_backslip", false,
      "Enable backslip loading: prescribe V=Vp below backslip_depth");
  params.addParam<Real>("Vp", 1e-9, "Plate rate for backslip loading (m/s)");
  params.addParam<Real>("backslip_depth", 18000.0,
      "Depth below which V=Vp is prescribed (m). Default: H+h=18km for BP2");
  params.addParam<unsigned int>("depth_direction", 1,
      "Coordinate direction for depth (0=x, 1=y, 2=z). Default: y");

  return params;
}

SEASSlipRateVarAAux::SEASSlipRateVarAAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _traction(coupledValue("traction")),
    _state_variable(coupledValue("state_variable")),
    _a_var(coupledValue("a_var")),
    _b(getParam<Real>("b")),
    _Dc(getParam<Real>("Dc")),
    _f0(getParam<Real>("f0")),
    _V0(getParam<Real>("V0")),
    _sigma_n(getParam<Real>("sigma_n")),
    _eta(getParam<Real>("eta")),
    _V_min(getParam<Real>("V_min")),
    _V_max(getParam<Real>("V_max")),
    _tau_pre(getParam<Real>("tau_pre")),
    _use_backslip(getParam<bool>("use_backslip")),
    _Vp(getParam<Real>("Vp")),
    _backslip_depth(getParam<Real>("backslip_depth")),
    _depth_dir(getParam<unsigned int>("depth_direction"))
{
}

Real
SEASSlipRateVarAAux::computeValue()
{
  // Backslip loading: if enabled, prescribe V = Vp below backslip_depth
  // This creates the correct stress loading pattern for half-space SEAS simulations
  // Note: depth coordinate may be negative (y < 0 for depth), so use abs()
  if (_use_backslip)
  {
    Real depth = std::fabs(_q_point[_qp](_depth_dir));
    if (depth > _backslip_depth)
      return _Vp;  // Prescribed plate rate in deep region
  }

  // Total traction = elastic traction + pre-stress
  Real tau = _traction[_qp] + _tau_pre;
  Real theta = _state_variable[_qp];
  Real a_local = _a_var[_qp];

  // Ensure positive values
  if (theta <= 0.0)
    theta = _Dc / _V_min;
  if (a_local <= 0.0)
    a_local = 0.025;  // Fallback to amax

  // Solve for slip rate from traction balance
  return solveSlipRate(tau, theta, a_local);
}

Real
SEASSlipRateVarAAux::solveSlipRate(Real tau, Real theta, Real a_local) const
{
  // For tension (sigma_n <= 0), slip rate is simply tau/eta
  if (_sigma_n <= 0.0)
    return std::max(std::fabs(tau) / _eta, _V_min);

  // Handle negative traction
  Real tau_abs = std::fabs(tau);

  // Define the residual function: R(V) = τ - σn * f(V, θ) - η * V
  auto residual = [&](double V) {
    Real f = frictionCoefficient(V, theta, a_local);
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
    return std::max(tau_abs / (a_local * _sigma_n / _V0 + _eta), _V_min);
  }

  // Ensure proper bracket ordering (R should go from + to -)
  if (R_min < R_max)
  {
    std::swap(V_min, V_max);
  }

  // Find the root using Brent's method
  try
  {
    Real V_solved = BrentRootFinder::zeroIn(V_min, V_max, residual);
    // Apply maximum slip rate limit for stability
    return std::min(V_solved, _V_max);
  }
  catch (const std::exception &)
  {
    // Fallback to linear approximation
    Real V_approx = std::max(tau_abs / (a_local * _sigma_n / _V0 + _eta), _V_min);
    return std::min(V_approx, _V_max);
  }
}

Real
SEASSlipRateVarAAux::frictionCoefficient(Real V, Real theta, Real a_local) const
{
  // Regularized friction coefficient (Tandem paper Eq. 7):
  // f(V, θ) = a * asinh[V/(2*V0) * exp((f0 + b*ln(V0*θ/Dc))/a)]

  // Ensure positive arguments for log
  V = std::max(V, _V_min);
  theta = std::max(theta, _Dc / _V0);

  Real arg = V / (2.0 * _V0) * std::exp((_f0 + _b * std::log(_V0 * theta / _Dc)) / a_local);
  return a_local * std::asinh(arg);
}
