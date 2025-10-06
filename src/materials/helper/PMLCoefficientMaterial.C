//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

/*
 * Material that computes PML damping coefficient
 */

#include "PMLCoefficientMaterial.h"

registerMooseObject("farmsApp", PMLCoefficientMaterial);

InputParameters
PMLCoefficientMaterial::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Computes spatially-varying PML damping coefficient "
                             "that increases from 0 at the physical domain edge "
                             "to d_max at the outer boundary.");

  // PML geometry parameters
  params.addRequiredParam<Real>("pml_xmin", "Inner x-coordinate of PML (left edge)");
  params.addRequiredParam<Real>("pml_xmax", "Inner x-coordinate of PML (right edge)");
  params.addRequiredParam<Real>("pml_ymin", "Inner y-coordinate of PML (bottom edge)");
  params.addRequiredParam<Real>("pml_ymax", "Inner y-coordinate of PML (top edge)");
  params.addRequiredParam<Real>("pml_thickness", "Thickness of PML layer");

  // Damping parameters
  params.addParam<Real>("d_max",
                        "Maximum damping coefficient. If not specified, computed as "
                        "3*ref_wave_speed/(2*pml_thickness) for ~1% reflection");
  params.addParam<Real>("exponent", 2.0, "Exponent for damping profile (typically 2-4)");
  params.addRequiredParam<Real>("ref_wave_speed",
                                "Reference wave speed for scaling (use P-wave speed)");

  return params;
}

PMLCoefficientMaterial::PMLCoefficientMaterial(const InputParameters & parameters)
  : Material(parameters),
    _pml_damping_coeff(declareProperty<Real>("pml_damping_coeff")),
    _pml_xmin(getParam<Real>("pml_xmin")),
    _pml_xmax(getParam<Real>("pml_xmax")),
    _pml_ymin(getParam<Real>("pml_ymin")),
    _pml_ymax(getParam<Real>("pml_ymax")),
    _pml_thickness(getParam<Real>("pml_thickness")),
    _d_max(isParamValid("d_max")
               ? getParam<Real>("d_max")
               : 3.0 * getParam<Real>("ref_wave_speed") / (2.0 * getParam<Real>("pml_thickness"))),
    _exponent(getParam<Real>("exponent")),
    _ref_wave_speed(getParam<Real>("ref_wave_speed"))
{
  // Validation
  if (_pml_thickness <= 0.0)
    mooseError("PML thickness must be positive");
  if (_d_max < 0.0)
    mooseError("Maximum damping coefficient must be non-negative");
  if (_exponent < 1.0)
    mooseError("Damping exponent should be >= 1.0");
}

void
PMLCoefficientMaterial::computeQpProperties()
{
  // Get current point coordinates
  const Point & p = _q_point[_qp];

  // Compute distance into PML
  Real r = computePMLDistance(p);

  // Compute damping coefficient using power law profile
  if (r <= 0.0)
  {
    // Inside physical domain - no damping
    _pml_damping_coeff[_qp] = 0.0;
  }
  else if (r >= _pml_thickness)
  {
    // Beyond PML outer boundary - maximum damping
    _pml_damping_coeff[_qp] = _d_max;
  }
  else
  {
    // Inside PML - gradual increase
    Real normalized_distance = r / _pml_thickness;
    _pml_damping_coeff[_qp] = _d_max * std::pow(normalized_distance, _exponent);
  }
}

Real
PMLCoefficientMaterial::computePMLDistance(const Point & p) const
{
  Real x = p(0);
  Real y = p(1);

  // Compute perpendicular distance from physical domain edge
  Real dist_x = 0.0;
  Real dist_y = 0.0;

  // X-direction distance
  if (x < _pml_xmin)
    dist_x = _pml_xmin - x;
  else if (x > _pml_xmax)
    dist_x = x - _pml_xmax;

  // Y-direction distance
  if (y < _pml_ymin)
    dist_y = _pml_ymin - y;
  else if (y > _pml_ymax)
    dist_y = y - _pml_ymax;

  // For corner regions, use maximum distance (conservative approach)
  // Alternative: use sqrt(dist_x^2 + dist_y^2) for diagonal distance
  return std::max(dist_x, dist_y);
}
