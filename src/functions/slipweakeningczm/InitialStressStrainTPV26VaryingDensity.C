#include "InitialStressStrainTPV26VaryingDensity.h"
#include <vector>
#include <algorithm>
#include <cmath>

registerMooseObject("farmsApp", InitialStressStrainTPV26VaryingDensity);

InputParameters
InitialStressStrainTPV26VaryingDensity::validParams()
{
  InputParameters params = Function::validParams();
  params.addClassDescription("Initial stress/strain (TPV26) using depth-varying density from TPV32 profile.");
  params.addRequiredParam<Real>("i", "index");
  params.addRequiredParam<Real>("j", "index");
  params.addRequiredParam<Real>("lambda_o", "initial lambda parameter for the CDBM model");
  params.addRequiredParam<Real>("shear_modulus_o", "initial shear modulus parameter for the CDBM model");
  params.addRequiredParam<Real>("fluid_density", "fluid density in kg/m^3");
  params.addRequiredParam<Real>("gravity", "gravity in m/s^2");
  params.addRequiredParam<Real>("bxx", "coefficient for sigmaxx");
  params.addRequiredParam<Real>("byy", "coefficient for sigmayy");
  params.addRequiredParam<Real>("bxy", "coefficient for sigmaxy");
  params.addParam<bool>("get_initial_stress", false, "flag to get initial stress");
  params.addParam<bool>("get_initial_strain", false, "flag to get initial strain");
  params.addParam<bool>("get_shear_overstress", false, "flag to get initial overstress");
  params.addParam<bool>("get_fluid_pressure", false, "flag to get fluid pressure");
  params.addParam<bool>("use_tapering", false, "use tapering to reduce deviatoric stress at shallow depth");
  params.addParam<Real>("tapering_depth_A", -1, "depth at which tapering starts");
  params.addParam<Real>("tapering_depth_B", -1, "depth at which tapering stops");
  params.addParam<bool>("use_overpressure", false, "use overpressure model");
  params.addParam<Real>("overpressure_depth_A", -1, "overpressure A depth");
  params.addParam<Real>("overpressure_depth_B", -1, "overpressure B depth");
  params.addParam<bool>("flip_sign", true, "If true, depth = -z (negative z is down)");
  params.addParam<Real>("depth_offset", 0.0, "Offset added to z before computing depth");
  return params;
}

InitialStressStrainTPV26VaryingDensity::InitialStressStrainTPV26VaryingDensity(const InputParameters & parameters)
  : Function(parameters),
    _i(getParam<Real>("i")),
    _j(getParam<Real>("j")),
    _lambda_o(getParam<Real>("lambda_o")),
    _shear_modulus_o(getParam<Real>("shear_modulus_o")),
    _fluid_density(getParam<Real>("fluid_density")),
    _gravity(getParam<Real>("gravity")),
    _bxx(getParam<Real>("bxx")),
    _byy(getParam<Real>("byy")),
    _bxy(getParam<Real>("bxy")),
    _get_initial_stress(getParam<bool>("get_initial_stress")),
    _get_initial_strain(getParam<bool>("get_initial_strain")),
    _get_shear_overstress(getParam<bool>("get_shear_overstress")),
    _get_fluid_pressure(getParam<bool>("get_fluid_pressure")),
    _use_tapering(getParam<bool>("use_tapering")),
    _tapering_depth_A(getParam<Real>("tapering_depth_A")),
    _tapering_depth_B(getParam<Real>("tapering_depth_B")),
    _use_overpressure(getParam<bool>("use_overpressure")),
    _overpressure_depth_A(getParam<Real>("overpressure_depth_A")),
    _overpressure_depth_B(getParam<Real>("overpressure_depth_B")),
    _flip_sign(getParam<bool>("flip_sign")),
    _depth_offset(getParam<Real>("depth_offset"))
{
  if (_get_initial_stress && _get_initial_strain)
    mooseError("Cannot get both initial stress and strain at the same time. Please choose one.");
  if (_get_fluid_pressure && (_get_initial_stress || _get_initial_strain))
    mooseError("When get_fluid_pressure is true, get_initial_stress and get_initial_strain must be false.");
  if (_use_tapering && (_tapering_depth_A < 0 || _tapering_depth_B < 0 || _tapering_depth_A >= _tapering_depth_B))
    mooseError("When use_tapering is true, tapering_depth_A and tapering_depth_B must be provided and A < B.");
  if (_use_overpressure && (_overpressure_depth_A < 0 || _overpressure_depth_B < 0 || _overpressure_depth_A >= _overpressure_depth_B))
    mooseError("When use_overpressure is true, overpressure_depth_A and overpressure_depth_B must be provided and A < B.");
}

Real InitialStressStrainTPV26VaryingDensity::tpv32_density_at_depth(Real d) const
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
InitialStressStrainTPV26VaryingDensity::value(Real /*t*/, const Point & p) const
{
  Real var = 0.0;

  const Real x = p(0);
  const Real y = p(1);
  const Real z = p(2);

  const Real depth = (_flip_sign ? -(z + _depth_offset) : (z + _depth_offset));
  const Real rock_density = tpv32_density_at_depth(std::max(depth, 0.0));

  // locals
  const Real lambda_o = _lambda_o;
  const Real shear_modulus_o = _shear_modulus_o;
  const Real fluid_density = _fluid_density;
  const Real gravity = _gravity;
  const Real bxx = _bxx;
  const Real byy = _byy;
  const Real bxy = _bxy;

  // Stresses
  Real sigmazz = 0.0, sigmaxx = 0.0, sigmayy = 0.0, sigmaxy = 0.0, sigmaxz = 0.0, sigmayz = 0.0;

  // Fluid pressure Pf
  Real Pf = 0.0;
  if (!_use_overpressure)
    Pf = fluid_density * gravity * std::abs(z);
  else
  {
    if (std::abs(z) <= _overpressure_depth_A)
      Pf = fluid_density * gravity * std::abs(z);
    else if (std::abs(z) > _overpressure_depth_A && std::abs(z) <= _overpressure_depth_B)
    {
      const Real Pf_A = fluid_density * gravity * _overpressure_depth_A;
      const Real delta_rho = rock_density - fluid_density;
      Pf = Pf_A + gravity * (fluid_density * (std::abs(z) - _overpressure_depth_A) +
            0.5 * delta_rho * std::pow((std::abs(z) - _overpressure_depth_A), 2) /
            (_overpressure_depth_B - _overpressure_depth_A));
    }
    else
    {
      const Real Pf_B = fluid_density * gravity * _overpressure_depth_B +
                        0.5 * gravity * (rock_density - fluid_density) *
                        (_overpressure_depth_B - _overpressure_depth_A);
      Pf = Pf_B + rock_density * gravity * (std::abs(z) - _overpressure_depth_B);
    }
  }

  // Tapering coefficient
  Real Omega = 1.0;
  if (_use_tapering && std::abs(z) < _tapering_depth_A)
    Omega = 1.0;
  else if (_use_tapering && std::abs(z) >= _tapering_depth_A && std::abs(z) <= _tapering_depth_B)
    Omega = (_tapering_depth_B - std::abs(z)) / (_tapering_depth_B - _tapering_depth_A);
  else if (_use_tapering && std::abs(z) > _tapering_depth_B)
    Omega = 0.0;

  // Vertical stress: use overburden integral with depth-varying density rho(z)
  // sigma_zz(z) = - g * \int_0^{depth} rho(s) ds, where depth = max( adjusted z, 0 )
  {
    const Real D = std::max(depth, 0.0);
    Real overburden = 0.0; // integral of rho(s) ds from 0 to D [kg/m^2]

    // TPV32 piecewise-linear density profile knots (depth in m, rho in kg/m^3)
    static const std::vector<Real> k_depth = {0, 500, 1000, 1600, 2400, 3600, 5000, 9000, 11000, 15000};
    static const std::vector<Real> k_rho   = {2200, 2450, 2550, 2600, 2600, 2620, 2650, 2720, 2750, 2900};

    if (D <= k_depth.front())
    {
      overburden = 0.0;
    }
    else
    {
      // Accumulate segment-wise using trapezoidal rule within each linear segment
      const std::size_t n = k_depth.size();
      for (std::size_t i = 0; i + 1 < n; ++i)
      {
        const Real z0 = k_depth[i];
        const Real z1 = k_depth[i + 1];
        if (D <= z0)
          break; // nothing more to add
        const Real upper = std::min(D, z1);
        if (upper > z0)
        {
          const Real rho0 = k_rho[i];
          const Real rho1 = k_rho[i + 1];
          const Real dz = upper - z0;
          // linear interpolation at 'upper'
          const Real rho_upper = rho0 + (rho1 - rho0) * (upper - z0) / (z1 - z0);
          overburden += 0.5 * (rho0 + rho_upper) * dz;
        }
        if (D <= z1)
          break; // finished within this segment
      }

      // If beyond last knot, extrapolate as constant density = k_rho.back()
      if (D > k_depth.back())
        overburden += k_rho.back() * (D - k_depth.back());
    }

    sigmazz = -gravity * overburden;
  }

  // Lateral and shear with tapering
  sigmaxx = Omega * (bxx * (sigmazz + Pf) - Pf) + (1.0 - Omega) * sigmazz;
  sigmayy = Omega * (byy * (sigmazz + Pf) - Pf) + (1.0 - Omega) * sigmazz;
  sigmaxy = Omega * (bxy * (sigmazz + Pf));

  //convert total stress to effective stress
  sigmaxx = sigmaxx + Pf;
  sigmayy = sigmayy + Pf;
  sigmazz = sigmazz + Pf;

  // Strain via inverse Hooke's law with constant lambda_o, shear_modulus_o (as in original)
  const Real sigma_mean = (sigmaxx + sigmayy + sigmazz);
  const Real first = (1.0 / (2.0 * shear_modulus_o));
  const Real second = (lambda_o / (2 * shear_modulus_o * (3 * lambda_o + 2 * shear_modulus_o)));
  const Real epsxx = first * sigmaxx - second * sigma_mean;
  const Real epsyy = first * sigmayy - second * sigma_mean;
  const Real epszz = first * sigmazz - second * sigma_mean;
  const Real epsxy = first * sigmaxy;
  const Real epsxz = first * sigmaxz;
  const Real epsyz = first * sigmayz;

  if (_get_initial_stress)
  {
    if (_i == 1 && _j == 1) var = sigmaxx;
    else if (_i == 2 && _j == 2) var = sigmayy;
    else if (_i == 3 && _j == 3) var = sigmazz;
    else if ((_i == 1 && _j == 2) || (_i == 2 && _j == 1)) var = sigmaxy;
    else if ((_i == 1 && _j == 3) || (_i == 3 && _j == 1)) var = sigmaxz;
    else if ((_i == 2 && _j == 3) || (_i == 3 && _j == 2)) var = sigmayz;
    else var = 0.0;
  }
  else if (_get_initial_strain)
  {
    if (_i == 1 && _j == 1) var = epsxx;
    else if (_i == 2 && _j == 2) var = epsyy;
    else if (_i == 3 && _j == 3) var = epszz;
    else if ((_i == 1 && _j == 2) || (_i == 2 && _j == 1)) var = epsxy;
    else if ((_i == 1 && _j == 3) || (_i == 3 && _j == 1)) var = epsxz;
    else if ((_i == 2 && _j == 3) || (_i == 3 && _j == 2)) var = epsyz;
    else var = 0.0;
  }
  else if (_get_fluid_pressure)
  {
    var = Pf;
  }

  return var;
}
