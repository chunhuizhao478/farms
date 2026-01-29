//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "RadiationDampingMaterial.h"
#include <cmath>

registerMooseObject("farmsApp", RadiationDampingMaterial);

InputParameters
RadiationDampingMaterial::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription(
      "Computes the radiation damping coefficient for quasi-dynamic SEAS simulations. "
      "η = μ / (2 * cs) where cs = sqrt(μ/ρ) is the shear wave speed.");

  params.addRequiredParam<Real>("shear_modulus", "Shear modulus μ (Pa)");
  params.addRequiredParam<Real>("density", "Mass density ρ (kg/m³)");

  params.addParam<MaterialPropertyName>("radiation_damping_name", "radiation_damping",
                                        "Name of the radiation damping property");
  params.addParam<MaterialPropertyName>("shear_wave_speed_name", "shear_wave_speed",
                                        "Name of the shear wave speed property");

  return params;
}

RadiationDampingMaterial::RadiationDampingMaterial(const InputParameters & parameters)
  : Material(parameters),
    _shear_modulus(getParam<Real>("shear_modulus")),
    _density(getParam<Real>("density")),
    _radiation_damping(
        declareProperty<Real>(getParam<MaterialPropertyName>("radiation_damping_name"))),
    _shear_wave_speed(
        declareProperty<Real>(getParam<MaterialPropertyName>("shear_wave_speed_name")))
{
  // Validate inputs
  if (_shear_modulus <= 0)
    paramError("shear_modulus", "Shear modulus must be positive");
  if (_density <= 0)
    paramError("density", "Density must be positive");
}

void
RadiationDampingMaterial::computeQpProperties()
{
  // Shear wave speed: cs = sqrt(μ/ρ)
  _shear_wave_speed[_qp] = std::sqrt(_shear_modulus / _density);

  // Radiation damping coefficient: η = μ / (2 * cs)
  _radiation_damping[_qp] = _shear_modulus / (2.0 * _shear_wave_speed[_qp]);
}
