//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ElkADPorousFlowDamagedBiotCoefficient.h"

#include <algorithm>
#include <cmath>

registerMooseObject("farmsApp", ElkADPorousFlowDamagedBiotCoefficient);

InputParameters
ElkADPorousFlowDamagedBiotCoefficient::validParams()
{
  InputParameters params = ADMaterial::validParams();
  params.addRequiredCoupledVar("phase_field", "Damage (phase-field) variable, d.");
  params.addRequiredRangeCheckedParam<Real>(
      "solid_bulk_compliance",
      "solid_bulk_compliance>0",
      "Reciprocal of the intact drained bulk modulus of the porous skeleton (1 / K_0).");
  params.addRequiredRangeCheckedParam<Real>(
      "grain_bulk_modulus", "grain_bulk_modulus>0", "Solid grain bulk modulus K_s");
  params.addParam<Real>(
      "minimum_degradation", 1e-8, "Lower bound applied to g(d) = (1 - d)^2 for stability.");
  params.addClassDescription("AD material that computes Biot coefficient alpha = 1 - K/K_s using damaged K = g(d)*K_0");
  return params;
}

ElkADPorousFlowDamagedBiotCoefficient::ElkADPorousFlowDamagedBiotCoefficient(
    const InputParameters & parameters)
  : ADMaterial(parameters),
    _damage(adCoupledValue("phase_field")),
    _bulk_modulus_intact(1.0 / getParam<Real>("solid_bulk_compliance")),
    _grain_bulk_modulus(getParam<Real>("grain_bulk_modulus")),
    _min_degradation(getParam<Real>("minimum_degradation")),
    _biot_coefficient_damaged(declareADProperty<Real>("biot_coefficient_damaged"))
{
}

void
ElkADPorousFlowDamagedBiotCoefficient::computeQpProperties()
{
  // Clamp damage to [0, 1]
  ADReal dmg = std::min(std::max(_damage[_qp], ADReal(0.0)), ADReal(1.0));

  // Degradation function g(d) = (1 - d)^2, with minimum for stability
  ADReal g = std::max(pow(1.0 - dmg, 2), ADReal(_min_degradation));

  // Damaged bulk modulus
  ADReal Kc = g * _bulk_modulus_intact;

  // Biot coefficient: alpha = 1 - K(d) / K_s
  ADReal alpha = 1.0 - Kc / _grain_bulk_modulus;

  // Clamp to [0, 1] for numerical safety
  _biot_coefficient_damaged[_qp] = std::min(std::max(alpha, ADReal(0.0)), ADReal(1.0));
}
