//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ElkPorousFlowDamagedBiotCoefficient.h"

#include <algorithm>
#include <cmath>

registerMooseObject("farmsApp", ElkPorousFlowDamagedBiotCoefficient);

InputParameters
ElkPorousFlowDamagedBiotCoefficient::validParams()
{
  InputParameters params = Material::validParams();
  params.addRequiredCoupledVar("phase_field", "Damage (phase-field) variable, c.");
  params.addRequiredRangeCheckedParam<Real>(
      "solid_bulk_compliance",
      "solid_bulk_compliance>0",
      "Reciprocal of the intact drained bulk modulus of the porous skeleton (1 / K_0).");
  params.addRequiredRangeCheckedParam<Real>(
      "grain_bulk_modulus", "grain_bulk_modulus>0", "Solid grain bulk modulus K_s");
  params.addParam<Real>(
      "minimum_degradation", 1e-8, "Lower bound applied to g(c) = (1 - c)^2 for stability.");
  params.addClassDescription("Computes Biot coefficient alpha = 1 - K/K_s using damaged K = 1/Cd");
  return params;
}

ElkPorousFlowDamagedBiotCoefficient::ElkPorousFlowDamagedBiotCoefficient(
    const InputParameters & parameters)
  : Material(parameters),
    _damage(coupledValue("phase_field")),
    _bulk_modulus_intact(1.0 / getParam<Real>("solid_bulk_compliance")),
    _grain_bulk_modulus(getParam<Real>("grain_bulk_modulus")),
    _min_degradation(getParam<Real>("minimum_degradation")),
    _biot_coefficient_damaged(declareProperty<Real>("biot_coefficient_damaged"))
{
}

void
ElkPorousFlowDamagedBiotCoefficient::computeQpProperties()
{
  const Real dmg = std::clamp(_damage[_qp], 0.0, 1.0);
  const Real g = std::max(std::pow(1.0 - dmg, 2), _min_degradation); //only valid for AT1 or AT2 model
  const Real Kc = g * _bulk_modulus_intact;
  const Real alpha = 1.0 - Kc / _grain_bulk_modulus;
  // note here the biot coefficient depends on the damage variable, only bulk modulus is modified Kc = g * _bulk_modulus_intact
  // clamp to [0,1] for numerical safety
  const Real alpha_clamped = std::clamp(alpha, 0.0, 1.0);
  _biot_coefficient_damaged[_qp] = alpha_clamped;
}
