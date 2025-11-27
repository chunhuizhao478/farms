#include "EffectiveBodyForceTPV26VaryingDensity.h"

registerMooseObject("farmsApp", EffectiveBodyForceTPV26VaryingDensity);

InputParameters
EffectiveBodyForceTPV26VaryingDensity::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Kernel to apply effective body force accounting for pore pressure gradient with depth-varying density from TPV32 profile: f_eff = rho(z) * g - dPf/dz. Extends EffectiveBodyForceTPV26 to support depth-varying density.");
  params.addRequiredParam<Real>("fluid_density", "fluid density in kg/m^3");
  params.addRequiredParam<Real>("gravity", "gravity in m/s^2");
  params.addParam<bool>("use_overpressure", false, "flag to use overpressure in the calculation, default is false");
  params.addParam<Real>("overpressure_depth_A", -1, "depth at which overpressure starts to transition");
  params.addParam<Real>("overpressure_depth_B", -1, "depth at which overpressure stops to transition");
  params.addParam<Real>("overpressure_rho_ref", 2900.0, "reference rock density for overpressure gradient (default: TPV32 final value = 2900 kg/m^3)");
  return params;
}

EffectiveBodyForceTPV26VaryingDensity::EffectiveBodyForceTPV26VaryingDensity(const InputParameters & parameters)
  : Kernel(parameters),
  _fluid_density(getParam<Real>("fluid_density")),
  _gravity(getParam<Real>("gravity")),
  _use_overpressure(getParam<bool>("use_overpressure")),
  _overpressure_depth_A(getParam<Real>("overpressure_depth_A")),
  _overpressure_depth_B(getParam<Real>("overpressure_depth_B")),
  _overpressure_rho_ref(getParam<Real>("overpressure_rho_ref"))
{
  //some checks for parameters
  if (_use_overpressure && (_overpressure_depth_A < 0 || _overpressure_depth_B < 0 || _overpressure_depth_A >= _overpressure_depth_B)) {
    mooseError("When use_overpressure is true, overpressure_depth_A and overpressure_depth_B must be provided and A must be less than B.");
  }
}

Real
EffectiveBodyForceTPV26VaryingDensity::tpv32_density_at_depth(Real d) const
{
  static const std::vector<Real> depth = {0, 500, 1000, 1600, 2400, 3600, 5000, 9000, 11000, 15000};
  static const std::vector<Real> rho   = {2200, 2450, 2550, 2600, 2600, 2620, 2650, 2720, 2750, 2900};
  if (d <= depth.front()) return rho.front();
  if (d >= depth.back())  return rho.back();
  for (std::size_t i = 0; i + 1 < depth.size(); ++i)
  {
    if (d >= depth[i] && d <= depth[i+1])
    {
      const Real t = (d - depth[i]) / (depth[i+1] - depth[i]);
      return rho[i] + t * (rho[i+1] - rho[i]);
    }
  }
  return rho.back();
}

Real
EffectiveBodyForceTPV26VaryingDensity::computeQpResidual()
{
  // Get the vertical coordinate (z-direction, compression negative)
  Real z_coord = std::abs(_q_point[_qp](2));

  // Get local rock density at this depth from TPV32 profile
  Real rock_density = tpv32_density_at_depth(z_coord);

  // Compute dPf/dz (pore pressure gradient)
  Real dPf_dz = 0.0;

  if (!_use_overpressure) {
    // Standard hydrostatic case
    // dPf/dz = rho_fluid * g
    dPf_dz = _fluid_density * _gravity;
  }
  else {
    // Overpressure model (linear transition, uses reference density for overpressure gradient)

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
      Real delta_rho = _overpressure_rho_ref - _fluid_density;
      dPf_dz = _gravity * (_fluid_density + delta_rho * (z_coord - _overpressure_depth_A) / (_overpressure_depth_B - _overpressure_depth_A));
    }
    // Region 3: Fully lithostatic (depth > B)
    // Pf = Pf_B + rho_ref * g * (z - B)
    // dPf/dz = rho_ref * g
    else if (z_coord > _overpressure_depth_B) {
      dPf_dz = _overpressure_rho_ref * _gravity;
    }
  }

  // Compute effective body force
  // f_eff = rho(z) * g - dPf/dz
  // Use local depth-varying density from TPV32 profile
  // Return negative value (compression negative convention)
  Real f_eff = -1.0 * (rock_density * _gravity - dPf_dz);

  // Return the residual contribution (body force * test function)
  return -f_eff * _test[_i][_qp]; // negative because moving to RHS
}
