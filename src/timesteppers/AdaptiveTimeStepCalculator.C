// In AdaptiveTimeStepCalculator.C
#include "AdaptiveTimeStepCalculator.h"

registerMooseObject("farmsApp", AdaptiveTimeStepCalculator);

InputParameters
AdaptiveTimeStepCalculator::validParams()
{
  InputParameters params = TimeStepper::validParams();
  params.addClassDescription("Calculates the adaptive time step for quasi-dynamic simulations based on stability criteria");
  
  // Required parameters
  params.addRequiredParam<Real>("cp", "P-wave velocity");
  params.addParam<Real>("a_prem", 0.1, "Proportionality index for stable time step of CFL");
  params.addRequiredParam<Real>("shear_modulus", "Shear modulus (mu)");
  params.addRequiredParam<Real>("permeability", "Permeability value (k_py)");
  params.addRequiredParam<Real>("biot_modulus", "Inverse of Biot modulus (c_o)");
  
  // Rate-and-state friction parameters
  params.addRequiredParam<Real>("a_o", "Direct effect parameter (a)");
  params.addRequiredParam<Real>("b_o", "State evolution parameter (b)");
  params.addRequiredParam<Real>("L", "Characteristic slip distance");
  params.addRequiredParam<Real>("normal_stress", "Normal stress on the fault");
  params.addParam<Real>("zeta_max", 0.5, "Maximum allowed value for zeta parameter");
  
  // Maximum slip rate and min element size postprocessor
  params.addRequiredParam<PostprocessorName>("dx_min", "Minimum element size in the mesh");
  params.addRequiredParam<PostprocessorName>("max_slip_rate", "Postprocessor that provides the maximum slip rate");
  
  return params;
}

AdaptiveTimeStepCalculator::AdaptiveTimeStepCalculator(const InputParameters & parameters)
  : TimeStepper(parameters),
    PostprocessorInterface(this),
    _cp(getParam<Real>("cp")),
    _a_prem(getParam<Real>("a_prem")),
    _shear_modulus(getParam<Real>("shear_modulus")),
    _permeability(getParam<Real>("permeability")),
    _biot_modulus(getParam<Real>("biot_modulus")),
    _dx_min(getPostprocessorValue("dx_min")),
    _a_o(getParam<Real>("a_o")),
    _b_o(getParam<Real>("b_o")),
    _L(getParam<Real>("L")),
    _normal_stress(getParam<Real>("normal_stress")),
    _zeta_max(getParam<Real>("zeta_max")),
    _max_slip_rate(getPostprocessorValue("max_slip_rate"))
{
}


void
AdaptiveTimeStepCalculator::init()
{
  // Calculate the fixed time increments at initialization
  _dt_exp = _a_prem * _dx_min / _cp;
  
  // Compute zeta parameter for quasi-dynamics (Lapusta et al. 2000)
  Real gama = 1.0;
  Real h_star = gama * _shear_modulus * _L / (_normal_stress * std::abs(_a_o - _b_o));
  Real kai = 0.25 * (_b_o - _a_o) / _a_o * std::pow((h_star / _dx_min - 1.0), 2) - h_star / _dx_min;
  
  if (kai > 0)
    _zeta = _a_o / ((_b_o - _a_o) * (h_star / _dx_min - 1.0));
  else
    _zeta = 1.0 - _dx_min / h_star;
  
  // Ensure zeta is within bounds
  _zeta = std::min(0.5 * _zeta, _zeta_max);
}

Real
AdaptiveTimeStepCalculator::computeInitialDT()
{
  // Start with the explicit time step
  return _dt_exp;
}

Real
AdaptiveTimeStepCalculator::computeDT()
{
  // Hydraulic diffusivity
  Real c = _permeability * _biot_modulus;
  
  // Diffusion time increment
  _dt_diff = (_dx_min) * (_dx_min) / c;
  
  // Use the provided maximum slip rate
  Real max_slip_rate = _max_slip_rate;
  
  // If slip rate is too small, use a reasonable default to avoid division by zero
  if (max_slip_rate < 1e-12)
    max_slip_rate = 1e-12;
  
  // Calculate adaptive time step
  Real dt_evo = _zeta * _L / max_slip_rate;
  Real dt_test = 2.0 * _dt;
  
  // Discretize to be a multiple of the explicit time step
  Real dt_inc = std::floor(dt_evo / _dt_exp);
  dt_evo = dt_inc * _dt_exp;
  
  // Choose the minimum time step based on various criteria
  Real dtmax = 1.0; // Maximum allowed time step
  Real dt_calculated = std::min({dt_test, dt_evo, 2.0 * _dt, _dt_diff, dtmax});
  dt_calculated = std::max(dt_calculated, _dt_exp); // Don't go below explicit time step

  return dt_calculated;
}