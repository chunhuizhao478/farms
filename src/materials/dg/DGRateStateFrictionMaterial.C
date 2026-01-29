//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DGRateStateFrictionMaterial.h"
#include "BrentRootFinder.h"
#include "FEProblemBase.h"
#include <cmath>

registerMooseObject("farmsApp", DGRateStateFrictionMaterial);

InputParameters
DGRateStateFrictionMaterial::validParams()
{
  InputParameters params = InterfaceMaterial::validParams();
  params.addClassDescription(
      "Rate-and-state friction material for DG SEAS simulations. "
      "Uses Tandem-style decoupled approach: solve for slip rate from "
      "traction balance using Brent's method.");

  params.addRequiredCoupledVar("displacement", "The displacement variable");

  // Rate-state parameters
  params.addRequiredParam<Real>("a", "Direct effect parameter a");
  params.addRequiredParam<Real>("b", "Evolution effect parameter b");
  params.addRequiredParam<Real>("Dc", "Critical slip distance Dc (m)");
  params.addParam<Real>("f0", 0.6, "Reference friction coefficient");
  params.addParam<Real>("V0", 1e-6, "Reference slip velocity (m/s)");
  params.addRequiredParam<Real>("sigma_n", "Normal stress, positive in compression (Pa)");

  // Material parameters
  params.addRequiredParam<Real>("shear_modulus", "Shear modulus μ (Pa)");
  params.addRequiredParam<Real>("density", "Mass density ρ (kg/m³)");

  // Initial conditions
  params.addParam<Real>("initial_slip_rate", 1e-9, "Initial slip rate V (m/s)");
  params.addParam<Real>("initial_state_variable", 0.0,
                        "Initial state variable θ. If 0, computed from steady state.");
  params.addParam<Real>("tau_pre", 0.0, "Pre-stress (Pa). If 0, computed from steady state.");

  // Output property names
  params.addParam<MaterialPropertyName>("fault_traction_name", "fault_traction",
                                        "Name of the fault traction property");
  params.addParam<MaterialPropertyName>("dtraction_dslip_name", "dtraction_dslip",
                                        "Name of the traction derivative property");

  return params;
}

DGRateStateFrictionMaterial::DGRateStateFrictionMaterial(const InputParameters & parameters)
  : InterfaceMaterial(parameters),
    _u(coupledValue("displacement")),
    _u_neighbor(coupledNeighborValue("displacement")),
    _grad_u(coupledGradient("displacement")),
    _grad_u_neighbor(coupledNeighborGradient("displacement")),
    _normals(_assembly.normals()),
    _a(getParam<Real>("a")),
    _b(getParam<Real>("b")),
    _Dc(getParam<Real>("Dc")),
    _f0(getParam<Real>("f0")),
    _V0(getParam<Real>("V0")),
    _sigma_n(getParam<Real>("sigma_n")),
    _shear_modulus(getParam<Real>("shear_modulus")),
    _density(getParam<Real>("density")),
    _initial_slip_rate(getParam<Real>("initial_slip_rate")),
    _fault_traction(
        declarePropertyByName<Real>(getParam<MaterialPropertyName>("fault_traction_name"))),
    _slip(declareProperty<Real>("slip")),
    _slip_rate(declareProperty<Real>("slip_rate")),
    _state_variable(declareProperty<Real>("state_variable")),
    _dtraction_dslip(
        declarePropertyByName<Real>(getParam<MaterialPropertyName>("dtraction_dslip_name"))),
    _friction_coefficient(declareProperty<Real>("friction_coefficient")),
    _state_variable_old(getMaterialPropertyOld<Real>("state_variable")),
    _slip_old(getMaterialPropertyOld<Real>("slip")),
    _slip_rate_old(getMaterialPropertyOld<Real>("slip_rate")),
    _tau_pre(getParam<Real>("tau_pre"))
{
  // Compute radiation damping coefficient: η = μ / (2 * cs)
  Real cs = std::sqrt(_shear_modulus / _density);
  _eta = _shear_modulus / (2.0 * cs);
}

void
DGRateStateFrictionMaterial::initQpStatefulProperties()
{
  // Initialize slip as jump in displacement
  _slip[_qp] = _u[_qp] - _u_neighbor[_qp];

  // Initialize slip rate from parameter
  _slip_rate[_qp] = _initial_slip_rate;

  // Initialize state variable
  Real theta_init = getParam<Real>("initial_state_variable");
  if (theta_init <= 0)
  {
    // Compute steady-state theta: at steady state, dθ/dt = 0
    // => 1 - V*θ/Dc = 0 => θ = Dc/V
    _state_variable[_qp] = _Dc / _initial_slip_rate;
  }
  else
  {
    _state_variable[_qp] = theta_init;
  }

  // Initialize friction coefficient and traction
  _friction_coefficient[_qp] = frictionCoefficient(_initial_slip_rate, _state_variable[_qp]);

  // Compute initial traction (tau = f * sigma_n + eta * V at steady state)
  if (_tau_pre > 0.0)
    _fault_traction[_qp] = _tau_pre;
  else
    _fault_traction[_qp] = _friction_coefficient[_qp] * _sigma_n + _eta * _initial_slip_rate;

  // Initialize derivative (used for Jacobian)
  Real dfric_dV = dFrictionDSlipRate(_initial_slip_rate, _state_variable[_qp]);
  _dtraction_dslip[_qp] = (dfric_dV * _sigma_n + _eta) * _shear_modulus / (_Dc / _initial_slip_rate);
}

void
DGRateStateFrictionMaterial::computeQpProperties()
{
  // Get current time step
  const Real dt = _fe_problem.dt();

  // Compute current slip from displacement jump
  _slip[_qp] = _u[_qp] - _u_neighbor[_qp];

  // Get old state variable
  Real theta_old = _state_variable_old[_qp];

  // Compute elastic traction from displacement gradients
  // In antiplane shear, tau = mu * du/dn where n is the normal direction
  Real tau_elastic = computeElasticTraction();

  // Add pre-stress contribution if specified
  if (_tau_pre > 0.0)
    tau_elastic += _tau_pre;

  // Solve for slip rate V from traction balance using Brent's method:
  // tau_elastic = sigma_n * f(V, theta) + eta * V
  Real V;
  if (dt > 0)
  {
    try
    {
      V = solveSlipRate(tau_elastic, theta_old);
    }
    catch (const std::exception & e)
    {
      // If root finding fails, fall back to previous slip rate
      V = _slip_rate_old[_qp];
    }
  }
  else
  {
    // At initial time, use prescribed initial slip rate
    V = _initial_slip_rate;
  }

  // Ensure slip rate is positive (or use absolute value for mode III)
  V = std::abs(V);
  if (V < 1e-20)
    V = 1e-20; // Regularization to avoid log(0)

  _slip_rate[_qp] = V;

  // Update state variable using aging law with backward Euler
  // dθ/dt = 1 - V*θ/Dc
  // θ_new = θ_old + dt * (1 - V*θ_new/Dc)
  // θ_new * (1 + V*dt/Dc) = θ_old + dt
  // θ_new = (θ_old + dt) / (1 + V*dt/Dc)
  if (dt > 0)
    _state_variable[_qp] = (theta_old + dt) / (1.0 + V * dt / _Dc);
  else
    _state_variable[_qp] = theta_old;

  // Compute friction coefficient with updated state
  _friction_coefficient[_qp] = frictionCoefficient(V, _state_variable[_qp]);

  // Compute fault traction (friction traction plus radiation damping)
  _fault_traction[_qp] = _friction_coefficient[_qp] * _sigma_n + _eta * V;

  // Compute derivative for Jacobian
  // The traction depends on V, which depends on the elastic traction tau_elastic
  // d(tau)/d(slip) = (d(tau)/dV) * (dV/d(tau_elastic)) * (d(tau_elastic)/d(slip))
  //
  // At equilibrium: tau_elastic = sigma_n * f(V, theta) + eta * V
  // dV/d(tau_elastic) = 1 / (sigma_n * df/dV + eta)
  //
  // d(tau_elastic)/d(slip) ≈ mu / h (stiffness from DG flux)
  //
  // For the interface kernel, we provide dtraction_dslip = d(tau)/d(slip_jump)
  Real dfric_dV = dFrictionDSlipRate(V, _state_variable[_qp]);
  Real denom = _sigma_n * dfric_dV + _eta;
  if (std::abs(denom) > 1e-20)
  {
    // Effective stiffness from DG discretization (approximate)
    Real h = 1.0; // Mesh size factor (would need element size for accuracy)
    Real dV_dtau = 1.0 / denom;
    Real dtau_dslip = _shear_modulus / h;
    _dtraction_dslip[_qp] = (dfric_dV * _sigma_n + _eta) * dV_dtau * dtau_dslip;
  }
  else
  {
    _dtraction_dslip[_qp] = _eta * _shear_modulus / _Dc;
  }
}

Real
DGRateStateFrictionMaterial::computeElasticTraction() const
{
  // In antiplane shear (mode III), the shear traction is:
  // tau = mu * du/dn
  // where n is the normal direction to the fault
  //
  // For the DG formulation, we compute the average of gradients on both sides:
  // tau = mu * 0.5 * (grad_u + grad_u_neighbor) . n
  //
  // In 2D antiplane (x,y) with fault normal in x-direction:
  // tau = mu * 0.5 * (du/dx + du_neighbor/dx) * n_x

  // Get the normal vector (points from element to neighbor)
  const RealVectorValue & n = _normals[_qp];

  // Average gradient normal component
  Real grad_u_n = _grad_u[_qp] * n;
  Real grad_u_neighbor_n = _grad_u_neighbor[_qp] * n;
  Real avg_grad_n = 0.5 * (grad_u_n + grad_u_neighbor_n);

  // Elastic traction: tau = mu * du/dn
  return _shear_modulus * avg_grad_n;
}

Real
DGRateStateFrictionMaterial::solveSlipRate(Real tau_elastic, Real theta) const
{
  // For tension (sigma_n <= 0), slip rate is simply tau/eta
  if (_sigma_n <= 0.0)
    return std::fabs(tau_elastic) / _eta;

  // Define the residual function: R(V) = τ_elastic - σn * f(V, θ) - η * V
  auto residual = [&](double V) {
    Real f = frictionCoefficient(V, theta);
    return tau_elastic - _sigma_n * f - _eta * V;
  };

  // Bracket for V: [V_min, V_max]
  // V_min: small positive to avoid log(0)
  // V_max: at maximum, all traction is viscous: tau = eta * V
  Real V_min = 1e-20;
  Real V_max = std::max(std::fabs(tau_elastic) / _eta, 1.0);

  // Check signs at boundaries to ensure we have a bracket
  Real R_min = residual(V_min);
  Real R_max = residual(V_max);

  // If same sign, expand the bracket or use fallback
  if (R_min * R_max > 0)
  {
    // Try larger upper bound
    V_max = std::max(V_max * 100.0, 10.0);
    R_max = residual(V_max);

    if (R_min * R_max > 0)
    {
      // If still same sign, use simple linear approximation
      // tau ≈ (a * sigma_n / V0) * V + eta * V  for small V
      // V ≈ tau / (a * sigma_n / V0 + eta)
      return std::fabs(tau_elastic) / (_a * _sigma_n / _V0 + _eta);
    }
  }

  // Ensure proper order
  if (R_min < 0 && R_max > 0)
  {
    std::swap(V_min, V_max);
    std::swap(R_min, R_max);
  }

  // Find the root using Brent's method
  return BrentRootFinder::zeroIn(V_min, V_max, residual);
}

Real
DGRateStateFrictionMaterial::frictionCoefficient(Real V, Real theta) const
{
  // Regularized friction coefficient (Tandem paper Eq. 7):
  // f(V,θ) = a * asinh[V/(2*V0) * exp((f0 + b*ln(V0*θ/Dc))/a)]
  //
  // Can be rewritten as:
  // f = a * asinh[V/(2*V0) * exp(f0/a) * (V0*θ/Dc)^(b/a)]

  Real arg = V / (2.0 * _V0) * std::exp((_f0 + _b * std::log(_V0 * theta / _Dc)) / _a);
  return _a * std::asinh(arg);
}

Real
DGRateStateFrictionMaterial::dFrictionDSlipRate(Real V, Real theta) const
{
  // Derivative of regularized friction w.r.t. slip rate V
  // f = a * asinh(g(V,θ))
  // where g = V/(2*V0) * exp((f0 + b*ln(V0*θ/Dc))/a)
  //
  // df/dV = a * (1/sqrt(1 + g^2)) * dg/dV
  // dg/dV = 1/(2*V0) * exp((f0 + b*ln(V0*θ/Dc))/a)
  //       = g / V
  //
  // df/dV = a * (g/V) / sqrt(1 + g^2)
  //       = a / (V * sqrt(1 + 1/g^2))

  Real g = V / (2.0 * _V0) * std::exp((_f0 + _b * std::log(_V0 * theta / _Dc)) / _a);
  return _a * g / (V * std::sqrt(1.0 + g * g));
}
