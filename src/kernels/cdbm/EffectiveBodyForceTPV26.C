#include "EffectiveBodyForceTPV26.h"

registerMooseObject("farmsApp", EffectiveBodyForceTPV26);

InputParameters
EffectiveBodyForceTPV26::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Kernel to apply effective body force accounting for pore pressure gradient: f_eff = rho * g - dPf/dz. Replaces standard BodyForce kernel when using overpressure models.");
  params.addRequiredParam<Real>("fluid_density", "fluid density in kg/m^3");
  params.addRequiredParam<Real>("rock_density", "rock density in kg/m^3");
  params.addRequiredParam<Real>("gravity", "gravity in m/s^2");
  params.addParam<bool>("use_overpressure", false, "flag to use overpressure in the calculation, default is false");
  params.addParam<Real>("overpressure_depth_A", -1, "depth at which overpressure starts to transition");
  params.addParam<Real>("overpressure_depth_B", -1, "depth at which overpressure stops to transition");
  params.addParam<bool>("overpressure_loweffective", false, "flag to use low effective stress overpressure (quadratic transition, lambda_pp scaling), default is false. Requires use_overpressure = true");
  params.addParam<Real>("lambda_pp", 0.9, "pore pressure ratio for low effective stress overpressure, default is 0.9");
  return params;
}

EffectiveBodyForceTPV26::EffectiveBodyForceTPV26(const InputParameters & parameters)
  : Kernel(parameters),
  _fluid_density(getParam<Real>("fluid_density")),
  _rock_density(getParam<Real>("rock_density")),
  _gravity(getParam<Real>("gravity")),
  _use_overpressure(getParam<bool>("use_overpressure")),
  _overpressure_depth_A(getParam<Real>("overpressure_depth_A")),
  _overpressure_depth_B(getParam<Real>("overpressure_depth_B")),
  _overpressure_loweffective(getParam<bool>("overpressure_loweffective")),
  _lambda_pp(getParam<Real>("lambda_pp"))
{
  //some checks for parameters
  if (_use_overpressure && (_overpressure_depth_A < 0 || _overpressure_depth_B < 0 || _overpressure_depth_A >= _overpressure_depth_B)) {
    mooseError("When use_overpressure is true, overpressure_depth_A and overpressure_depth_B must be provided and A must be less than B.");
  }
  if (_overpressure_loweffective && !_use_overpressure) {
    mooseError("When overpressure_loweffective is true, use_overpressure must also be true.");
  }
  if (_overpressure_loweffective && (_lambda_pp <= 0.0 || _lambda_pp > 1.0)) {
    mooseError("When overpressure_loweffective is true, lambda_pp must be in the range (0, 1]. Typical values: 0.9 (10% effective stress), 0.95 (5%), 0.98 (2%).");
  }
}

Real
EffectiveBodyForceTPV26::computeQpResidual()
{
  // Get the vertical coordinate (z-direction, compression negative)
  Real z_coord = std::abs(_q_point[_qp](2));

  // Compute dPf/dz (pore pressure gradient)
  Real dPf_dz = 0.0;

  if (!_use_overpressure) {
    // Standard hydrostatic case
    // dPf/dz = rho_fluid * g
    dPf_dz = _fluid_density * _gravity;
  }
  else if (_overpressure_loweffective) {
    // Low effective stress overpressure model (quadratic transition, lambda_pp scaling)

    // Region 1: Hydrostatic (depth <= A)
    // Pf = rho_fluid * g * z
    // dPf/dz = rho_fluid * g
    if (z_coord <= _overpressure_depth_A) {
      dPf_dz = _fluid_density * _gravity;
    }
    // Region 2: Quadratic transition (A < depth <= B)
    // Pf = Pf_A + (Pf_B_target - Pf_A) * s^2, where s = (z - A) / (B - A)
    // dPf/dz = 2 * (Pf_B_target - Pf_A) / (B - A) * s
    else if (z_coord > _overpressure_depth_A && z_coord <= _overpressure_depth_B) {
      Real Pf_A = _fluid_density * _gravity * _overpressure_depth_A;
      Real Pf_B_target = _lambda_pp * _rock_density * _gravity * _overpressure_depth_B;
      Real s = (z_coord - _overpressure_depth_A) / (_overpressure_depth_B - _overpressure_depth_A);
      dPf_dz = 2.0 * (Pf_B_target - Pf_A) / (_overpressure_depth_B - _overpressure_depth_A) * s;
    }
    // Region 3: Scaled lithostatic (depth > B)
    // Pf = lambda_pp * rho_rock * g * z
    // dPf/dz = lambda_pp * rho_rock * g
    else if (z_coord > _overpressure_depth_B) {
      dPf_dz = _lambda_pp * _rock_density * _gravity;
    }
  }
  else {
    // Standard overpressure model (linear transition, full lithostatic)

    // Region 1: Hydrostatic (depth <= A)
    // Pf = rho_fluid * g * z
    // dPf/dz = rho_fluid * g
    if (z_coord <= _overpressure_depth_A) {
      dPf_dz = _fluid_density * _gravity;
    }
    // Region 2: Linear-gradient transition (A < depth <= B)
    // Pf = Pf_A + g * (rho_fluid * (z - A) + 0.5 * delta_rho * (z - A)^2 / (B - A))
    // dPf/dz = g * (rho_fluid + delta_rho * (z - A) / (B - A))
    else if (z_coord > _overpressure_depth_A && z_coord <= _overpressure_depth_B) {
      Real delta_rho = _rock_density - _fluid_density;
      dPf_dz = _gravity * (_fluid_density + delta_rho * (z_coord - _overpressure_depth_A) / (_overpressure_depth_B - _overpressure_depth_A));
    }
    // Region 3: Fully lithostatic (depth > B)
    // Pf = Pf_B + rho_rock * g * (z - B)
    // dPf/dz = rho_rock * g
    else if (z_coord > _overpressure_depth_B) {
      dPf_dz = _rock_density * _gravity;
    }
  }

  // Compute effective body force
  // f_eff = rho_rock * g - dPf/dz
  // Return negative value (compression negative convention)
  Real f_eff = -1.0 * (_rock_density * _gravity - dPf_dz);

  // Return the residual contribution (body force * test function)
  return -f_eff * _test[_i][_qp]; // negative because moving to RHS
}
