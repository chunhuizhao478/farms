//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ElkPorousFlowDamagedBiotModulus.h"

#include "MooseError.h"
#include "MooseException.h"

#include <algorithm>
#include <cmath>
#include <string>

registerMooseObject("farmsApp", ElkPorousFlowDamagedBiotModulus);

InputParameters
ElkPorousFlowDamagedBiotModulus::validParams()
{
  InputParameters params = PorousFlowMaterialVectorBase::validParams();
  params.addRangeCheckedParam<Real>(
      "biot_coefficient", 0.4, "biot_coefficient>=0 & biot_coefficient<=1", "Biot coefficient (constant, ignored if use_damaged_biot=true)");
  params.addRangeCheckedParam<Real>(
      "fluid_bulk_modulus", 2.0E9, "fluid_bulk_modulus>0", "Fluid bulk modulus");
  params.addRequiredRangeCheckedParam<Real>(
      "grain_bulk_modulus", "grain_bulk_modulus>0", "Solid grain bulk modulus (K_s)");
  params.addParam<bool>("use_damaged_biot", false, "Use biot coefficient from material property 'biot_coefficient_damaged'.");
  params.addParam<bool>("use_damaged_porosity", false, "Use porosity from material property 'PorousFlow_porosity_*_damaged' instead of constant value.");
  params.addRangeCheckedParam<Real>(
      "porosity",
      0.008,
      "porosity>=0 & porosity<=1",
      "Constant porosity used when use_damaged_porosity=false.");
  params.addPrivateParam<std::string>("pf_material_type", "biot_modulus");
  params.addClassDescription("Computes the damage-dependent Biot modulus using the porosity and "
                             "Biot coefficient relationships from the coupled phase-field model.");
  return params;
}

ElkPorousFlowDamagedBiotModulus::ElkPorousFlowDamagedBiotModulus(const InputParameters & parameters)
  : PorousFlowMaterialVectorBase(parameters),
    _biot_coefficient_const(getParam<Real>("biot_coefficient")),
    _use_damaged_biot(getParam<bool>("use_damaged_biot")),
    _biot_coefficient_damaged_matprop(nullptr),
    _fluid_bulk_modulus(getParam<Real>("fluid_bulk_modulus")),
    _grain_bulk_modulus(getParam<Real>("grain_bulk_modulus")),
    _use_damaged_porosity(getParam<bool>("use_damaged_porosity")),
    _porosity_damaged_matprop(nullptr),
    _porosity_const(getParam<Real>("porosity")),
    _biot_modulus(_nodal_material ? declareProperty<Real>("PorousFlow_constant_biot_modulus_nodal")
                                  : declareProperty<Real>("PorousFlow_constant_biot_modulus_qp"))
{
  if (_use_damaged_biot)
  {
    try
    {
      _biot_coefficient_damaged_matprop = &getMaterialProperty<Real>("biot_coefficient_damaged");
    }
    catch (const MooseException &)
    {
      mooseError("Requested damaged Biot coefficient but material property 'biot_coefficient_damaged' "
                 "was not found for ",
                 name(),
                 ".");
    }
  }

  if (_use_damaged_porosity)
  {
    try
    {
      _porosity_damaged_matprop =
          _nodal_material ? &getMaterialProperty<Real>("PorousFlow_porosity_nodal_damaged")
                          : &getMaterialProperty<Real>("PorousFlow_porosity_qp_damaged");
    }
    catch (const MooseException &)
    {
      mooseError("Requested damaged porosity but the property 'PorousFlow_porosity_*_damaged' "
                 "was not found for ",
                 name(),
                 ".");
    }
  }
  else
  {
    // Using constant porosity - validate range
    if (_porosity_const < 0.0 || _porosity_const > 1.0)
      mooseError("Parameter 'porosity' must be within [0,1] when use_damaged_porosity=false in ",
                 name(),
                 ".");
  }
}

void
ElkPorousFlowDamagedBiotModulus::initQpStatefulProperties()
{
  computeQpProperties();
}

void
ElkPorousFlowDamagedBiotModulus::computeQpProperties()
{
  //get biot coefficient
  const Real alpha = _use_damaged_biot ? (*_biot_coefficient_damaged_matprop)[_qp]
                                       : _biot_coefficient_const;

  //get porosity
  const Real phi = _use_damaged_porosity ? (*_porosity_damaged_matprop)[_qp]
                                         : _porosity_const;

  const Real denom =
      phi / _fluid_bulk_modulus + (alpha - phi) / _grain_bulk_modulus;
  const Real safe_denom = std::max(denom, _denominator_floor);
  _biot_modulus[_qp] = 1.0 / safe_denom;
}
