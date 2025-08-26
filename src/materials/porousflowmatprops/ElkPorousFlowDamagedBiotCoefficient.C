//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ElkPorousFlowDamagedBiotCoefficient.h"

registerMooseObject("farmsApp", ElkPorousFlowDamagedBiotCoefficient);

InputParameters
ElkPorousFlowDamagedBiotCoefficient::validParams()
{
  InputParameters params = Material::validParams();
  params.addRequiredRangeCheckedParam<Real>(
      "grain_bulk_modulus", "grain_bulk_modulus>0", "Solid grain bulk modulus K_s");
  params.addClassDescription("Computes Biot coefficient alpha = 1 - K/K_s using damaged K = 1/Cd");
  return params;
}

ElkPorousFlowDamagedBiotCoefficient::ElkPorousFlowDamagedBiotCoefficient(
    const InputParameters & parameters)
  : Material(parameters),
    _grain_bulk_modulus(getParam<Real>("grain_bulk_modulus")),
    _solid_bulk_compliance_damaged(getMaterialProperty<Real>("solid_bulk_compliance_damaged")),
    _biot_coefficient(declareProperty<Real>("biot_coefficient"))
{
}

void
ElkPorousFlowDamagedBiotCoefficient::computeQpProperties()
{
  const Real K = 1.0 / std::max(_solid_bulk_compliance_damaged[_qp], 1e-24);
  _biot_coefficient[_qp] = 1.0 - K / _grain_bulk_modulus;
  // clamp to [0,1] for numerical safety
  if (_biot_coefficient[_qp] < 0.0)
    _biot_coefficient[_qp] = 0.0;
  else if (_biot_coefficient[_qp] > 1.0)
    _biot_coefficient[_qp] = 1.0;
}
