//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SochackiSpongeMaterial.h"

#include <algorithm>
#include <cmath>

registerMooseObject("farmsApp", SochackiSpongeMaterial);

InputParameters
SochackiSpongeMaterial::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription(
      "Computes the Sochacki sponge attenuation coefficient A(x, y) used to absorb "
      "waves in the outer sponge layer. The coefficient ramps from zero at the "
      "interface between the physical domain and the sponge to S_max at the outer boundary.");

  params.addRequiredParam<Real>("inner_xmin",
                                "Minimum x-coordinate of the physical domain (start of sponge).");
  params.addRequiredParam<Real>("inner_xmax",
                                "Maximum x-coordinate of the physical domain (start of sponge).");
  params.addRequiredParam<Real>("inner_ymin",
                                "Minimum y-coordinate of the physical domain (start of sponge).");
  params.addRequiredParam<Real>("inner_ymax",
                                "Maximum y-coordinate of the physical domain (start of sponge).");

  params.addRequiredParam<Real>("outer_xmin",
                                "Minimum x-coordinate of the full mesh (outer sponge edge).");
  params.addRequiredParam<Real>("outer_xmax",
                                "Maximum x-coordinate of the full mesh (outer sponge edge).");
  params.addRequiredParam<Real>("outer_ymin",
                                "Minimum y-coordinate of the full mesh (outer sponge edge).");
  params.addRequiredParam<Real>("outer_ymax",
                                "Maximum y-coordinate of the full mesh (outer sponge edge).");

  params.addRequiredRangeCheckedParam<Real>(
      "s_max", "s_max>=0.0", "Maximum attenuation coefficient reached at the outer sponge edge.");

  MooseEnum profiles("linear exponent cubic exponential gaussian", "linear");
  params.addParam<MooseEnum>(
      "profile", profiles, "Shape of the damping ramp across the sponge thickness.");

  params.addParam<Real>(
      "exponent_power",
      2.0,
      "Power b used by the 'exponent' profile: A = S_max * xi^b, with xi in [0, 1].");
  params.addParam<Real>(
      "exp_rate",
      3.0,
      "Rate parameter used by the 'exponential' profile: exp_rate controls how fast the "
      "attenuation grows. A = S_max * (exp(exp_rate * xi) - 1)/(exp(exp_rate) - 1).");
  params.addParam<Real>(
      "gaussian_rate",
      3.0,
      "Rate parameter used by the 'gaussian' profile: A = S_max * (1 - exp(-gaussian_rate * xi^2))/"
      "(1 - exp(-gaussian_rate)).");

  return params;
}

SochackiSpongeMaterial::SochackiSpongeMaterial(const InputParameters & parameters)
  : Material(parameters),
    _sponge_coeff(declareProperty<Real>("sochacki_damping")),
    _inner_xmin(getParam<Real>("inner_xmin")),
    _inner_xmax(getParam<Real>("inner_xmax")),
    _inner_ymin(getParam<Real>("inner_ymin")),
    _inner_ymax(getParam<Real>("inner_ymax")),
    _outer_xmin(getParam<Real>("outer_xmin")),
    _outer_xmax(getParam<Real>("outer_xmax")),
    _outer_ymin(getParam<Real>("outer_ymin")),
    _outer_ymax(getParam<Real>("outer_ymax")),
    _s_max(getParam<Real>("s_max")),
    _exponent_power(getParam<Real>("exponent_power")),
    _exp_rate(getParam<Real>("exp_rate")),
    _gaussian_rate(getParam<Real>("gaussian_rate")),
    _profile(getParam<MooseEnum>("profile")),
    _left_thickness(_inner_xmin - _outer_xmin),
    _right_thickness(_outer_xmax - _inner_xmax),
    _bottom_thickness(_inner_ymin - _outer_ymin),
    _top_thickness(_outer_ymax - _inner_ymax),
    _exp_norm(_exp_rate > 0.0 ? std::exp(_exp_rate) - 1.0 : 1.0),
    _gaussian_norm(_gaussian_rate > 0.0 ? 1.0 - std::exp(-_gaussian_rate) : 1.0)
{
  if (_inner_xmax <= _inner_xmin || _inner_ymax <= _inner_ymin)
    mooseError("Inner bounds in SochackiSpongeMaterial must define a valid box.");

  if (_outer_xmax <= _outer_xmin || _outer_ymax <= _outer_ymin)
    mooseError("Outer bounds in SochackiSpongeMaterial must define a valid box.");

  if (_inner_xmin < _outer_xmin || _inner_xmax > _outer_xmax || _inner_ymin < _outer_ymin ||
      _inner_ymax > _outer_ymax)
    mooseError("Inner bounds must lie inside the outer bounds for Sochacki sponge.");

  if (_left_thickness < 0.0 || _right_thickness < 0.0 || _bottom_thickness < 0.0 ||
      _top_thickness < 0.0)
    mooseError("Computed sponge thickness is negative; check the supplied bounds.");

  if (_exponent_power < 1.0)
    mooseWarning("Exponent power < 1.0 may lead to non-monotonic attenuation.");
}

void
SochackiSpongeMaterial::computeQpProperties()
{
  const Real xi = computeNormalizedDistance(_q_point[_qp]);

  if (xi <= 0.0)
  {
    _sponge_coeff[_qp] = 0.0;
    return;
  }

  Real value = 0.0;

  const std::string & profile = _profile;
  if (profile == "linear")
    value = _s_max * xi;
  else if (profile == "exponent")
    value = _s_max * std::pow(xi, _exponent_power);
  else if (profile == "cubic")
    value = _s_max * xi * xi * xi;
  else if (profile == "exponential")
  {
    if (_exp_rate <= 0.0 || std::abs(_exp_norm) < 1e-12)
      value = _s_max * xi;
    else
      value = _s_max * (std::exp(_exp_rate * xi) - 1.0) / _exp_norm;
  }
  else if (profile == "gaussian")
  {
    if (_gaussian_rate <= 0.0 || std::abs(_gaussian_norm) < 1e-12)
      value = _s_max * xi;
    else
    {
      const Real num = 1.0 - std::exp(-_gaussian_rate * xi * xi);
      value = _s_max * num / _gaussian_norm;
    }
  }
  else
    mooseError("Unsupported Sochacki damping profile: ", profile);

  _sponge_coeff[_qp] = std::min(value, _s_max);
}

Real
SochackiSpongeMaterial::computeNormalizedDistance(const Point & p) const
{
  Real nx = 0.0;
  if (p(0) < _inner_xmin && _left_thickness > 0.0)
    nx = (_inner_xmin - p(0)) / _left_thickness;
  else if (p(0) > _inner_xmax && _right_thickness > 0.0)
    nx = (p(0) - _inner_xmax) / _right_thickness;

  Real ny = 0.0;
  if (p(1) < _inner_ymin && _bottom_thickness > 0.0)
    ny = (_inner_ymin - p(1)) / _bottom_thickness;
  else if (p(1) > _inner_ymax && _top_thickness > 0.0)
    ny = (p(1) - _inner_ymax) / _top_thickness;

  if (nx <= 0.0 && ny <= 0.0)
    return 0.0;

  nx = std::clamp(nx, 0.0, 1.0);
  ny = std::clamp(ny, 0.0, 1.0);

  if (nx > 0.0 && ny > 0.0)
    return std::clamp(std::sqrt(nx * ny), 0.0, 1.0);

  return std::max(nx, ny);
}
