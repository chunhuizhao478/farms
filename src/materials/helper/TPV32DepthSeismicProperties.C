//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "TPV32DepthSeismicProperties.h"

registerMooseObject("farmsApp", TPV32DepthSeismicProperties);

InputParameters
TPV32DepthSeismicProperties::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Depth-varying Vs, Vp, rho according to TPV32; outputs lambda_input and shear_modulus_input.");
  params.addParam<std::string>("depth_axis", "y", "Coordinate axis used as depth: x, y, or z");
  params.addParam<Real>("depth_offset", 0.0, "Offset added to coordinate before evaluating depth (meters)");
  params.addParam<bool>("clamp_nonpositive_depth", true, "Use 0 m properties for depth <= 0");
  params.addParam<bool>("flip_sign", false, "If true, depth = -(coord+offset) to support negative-down meshes");
  // Configurable output property names (backward-compatible defaults)
  params.addParam<MaterialPropertyName>("vs_property_name", "Vs", "Output property name for Vs (m/s)");
  params.addParam<MaterialPropertyName>("vp_property_name", "Vp", "Output property name for Vp (m/s)");
  params.addParam<MaterialPropertyName>("density_property_name", "density_input", "Output property name for density (kg/m^3)");
  params.addParam<MaterialPropertyName>("shear_modulus_property_name", "shear_modulus_input", "Output property name for shear modulus mu (Pa)");
  params.addParam<MaterialPropertyName>("lambda_property_name", "lambda_input", "Output property name for Lam\xC3\xA9's first parameter lambda (Pa)");
  return params;
}

TPV32DepthSeismicProperties::TPV32DepthSeismicProperties(const InputParameters & parameters)
  : Material(parameters),
    _axis([&]() {
      const auto ax = getParam<std::string>("depth_axis");
      if (ax == "x") return 0u; if (ax == "y") return 1u; if (ax == "z") return 2u;
      mooseError("depth_axis must be one of x,y,z");
    }()),
    _offset(getParam<Real>("depth_offset")),
  _clamp_nonpositive(getParam<bool>("clamp_nonpositive_depth")),
  _flip_sign(getParam<bool>("flip_sign")),
  _Vs(declarePropertyByName<Real>(getParam<MaterialPropertyName>("vs_property_name"))),
  _Vp(declarePropertyByName<Real>(getParam<MaterialPropertyName>("vp_property_name"))),
  _rho(declarePropertyByName<Real>(getParam<MaterialPropertyName>("density_property_name"))),
  _mu(declarePropertyByName<Real>(getParam<MaterialPropertyName>("shear_modulus_property_name"))),
  _lambda(declarePropertyByName<Real>(getParam<MaterialPropertyName>("lambda_property_name")))
{
}

Real TPV32DepthSeismicProperties::interp(const std::vector<Real> & xs,
                                         const std::vector<Real> & ys,
                                         Real x) const
{
  mooseAssert(xs.size() == ys.size(), "interp: size mismatch");
  if (x <= xs.front())
    return ys.front();
  if (x >= xs.back())
    return ys.back();
  // find interval [i,i+1] such that xs[i] <= x <= xs[i+1]
  for (std::size_t i = 0; i + 1 < xs.size(); ++i)
  {
    if (x >= xs[i] && x <= xs[i + 1])
    {
      const Real x0 = xs[i], x1 = xs[i + 1];
      const Real y0 = ys[i], y1 = ys[i + 1];
      const Real t = (x - x0) / (x1 - x0);
      return y0 + t * (y1 - y0);
    }
  }
  return ys.back();
}

void TPV32DepthSeismicProperties::computeQpProperties()
{
  // TPV32 table (depth in meters)
  static const std::vector<Real> depth = {0, 500, 1000, 1600, 2400, 3600, 5000, 9000, 11000, 15000};
  static const std::vector<Real> Vp = {2200, 3000, 3600, 4400, 4800, 5250, 5500, 5750, 6100, 6300};
  static const std::vector<Real> Vs = {1050, 1400, 1950, 2500, 2800, 3100, 3250, 3450, 3600, 3700};
  static const std::vector<Real> rho = {2200, 2450, 2550, 2600, 2600, 2620, 2650, 2720, 2750, 2900};

  // Depth calculation
  const Point & X = _q_point[_qp];
  Real d = X(_axis) + _offset;
  if (_flip_sign)
    d = -d;
  if (_clamp_nonpositive && d <= 0.0)
    d = 0.0;

  // Interpolate
  const Real vp = interp(depth, Vp, d);
  const Real vs = interp(depth, Vs, d);
  const Real r = interp(depth, rho, d);

  _Vp[_qp] = vp;
  _Vs[_qp] = vs;
  _rho[_qp] = r;

  // Derived elastic constants (SI units expected: m/s, kg/m^3)
  const Real mu = r * vs * vs;
  const Real lambda = std::max(r * vp * vp - 2.0 * mu, 0.0);

  _mu[_qp] = mu;
  _lambda[_qp] = lambda;
}
