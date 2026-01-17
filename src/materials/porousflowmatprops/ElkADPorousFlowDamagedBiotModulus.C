//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ElkADPorousFlowDamagedBiotModulus.h"

#include "MooseError.h"

#include <algorithm>
#include <cmath>

registerMooseObject("farmsApp", ElkADPorousFlowDamagedBiotModulus);

InputParameters
ElkADPorousFlowDamagedBiotModulus::validParams()
{
  InputParameters params = ADMaterial::validParams();
  params.addRangeCheckedParam<Real>(
      "biot_coefficient", 1.0, "biot_coefficient>=0 & biot_coefficient<=1",
      "Biot coefficient (constant, ignored if use_damaged_biot=true)");
  params.addRangeCheckedParam<Real>(
      "fluid_bulk_modulus", 2.0E9, "fluid_bulk_modulus>0", "Fluid bulk modulus");
  params.addRequiredRangeCheckedParam<Real>(
      "grain_bulk_modulus", "grain_bulk_modulus>0", "Solid grain bulk modulus (K_s)");
  params.addParam<bool>("use_damaged_biot", false,
      "Use biot coefficient from AD material property 'biot_coefficient_damaged'.");
  params.addParam<bool>("use_damaged_porosity", false,
      "Use porosity from AD material property 'PorousFlow_porosity_qp_damaged' instead of constant value.");
  params.addRangeCheckedParam<Real>(
      "porosity",
      0.008,
      "porosity>=0 & porosity<=1",
      "Constant porosity used when use_damaged_porosity=false.");
  params.addClassDescription("AD material that computes damage-dependent Biot modulus: 1/M = phi/Kf + (alpha-phi)/Ks");
  return params;
}

ElkADPorousFlowDamagedBiotModulus::ElkADPorousFlowDamagedBiotModulus(const InputParameters & parameters)
  : ADMaterial(parameters),
    _use_damaged_biot(getParam<bool>("use_damaged_biot")),
    _biot_coefficient_const(getParam<Real>("biot_coefficient")),
    _biot_coefficient_damaged_matprop(nullptr),
    _fluid_bulk_modulus(getParam<Real>("fluid_bulk_modulus")),
    _grain_bulk_modulus(getParam<Real>("grain_bulk_modulus")),
    _use_damaged_porosity(getParam<bool>("use_damaged_porosity")),
    _porosity_damaged_matprop(nullptr),
    _porosity_const(getParam<Real>("porosity")),
    _biot_modulus(declareADProperty<Real>(isNodal() ? "PorousFlow_constant_biot_modulus_nodal"
                                                    : "PorousFlow_constant_biot_modulus_qp"))
{
  if (_use_damaged_biot)
    _biot_coefficient_damaged_matprop = &getADMaterialProperty<Real>("biot_coefficient_damaged");

  if (_use_damaged_porosity)
    _porosity_damaged_matprop = &getADMaterialProperty<Real>(isNodal() ? "PorousFlow_porosity_nodal_damaged"
                                                                       : "PorousFlow_porosity_qp_damaged");
}

void
ElkADPorousFlowDamagedBiotModulus::computeQpProperties()
{
  // Get Biot coefficient (either from damaged property or constant)
  ADReal alpha = _use_damaged_biot ? (*_biot_coefficient_damaged_matprop)[_qp]
                                   : ADReal(_biot_coefficient_const);

  // Get porosity (either from damaged property or constant)
  ADReal phi = _use_damaged_porosity ? (*_porosity_damaged_matprop)[_qp]
                                     : ADReal(_porosity_const);

  // Compute 1/M = phi/Kf + (alpha - phi)/Ks
  ADReal denom = phi / _fluid_bulk_modulus + (alpha - phi) / _grain_bulk_modulus;

  // Apply floor for numerical stability
  ADReal safe_denom = std::max(denom, ADReal(_denominator_floor));

  // Biot modulus M = 1 / denom
  _biot_modulus[_qp] = 1.0 / safe_denom;
}
