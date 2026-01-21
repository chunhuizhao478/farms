//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ElkADPoroDynamicPhaseFieldMaterials.h"

/**
 * ADMaterial used in three-field poro-dynamics phase-field simulations.
 */
registerMooseObject("farmsApp", ElkADPoroDynamicPhaseFieldMaterials);

InputParameters
ElkADPoroDynamicPhaseFieldMaterials::validParams()
{
  InputParameters params = ADMaterial::validParams();
  params.addClassDescription("ADMaterial for three-field poro-dynamics with phase-field damage coupling");

  // Required parameters
  params.addRequiredParam<Real>("rhos_value", "Solid density value (kg/m^3)");
  params.addRequiredParam<Real>("rhof_value", "Fluid density value (kg/m^3)");
  params.addRequiredParam<Real>("tortosity_value", "Tortosity value");
  params.addRequiredParam<Real>("viscosity_value", "Fluid dynamic viscosity (Pa·s)");
  params.addRequiredParam<Real>("grain_bulk_modulus", "Solid grain bulk modulus K_s (Pa)");
  params.addRequiredParam<Real>("fluid_bulk_modulus", "Fluid bulk modulus K_f (Pa)");

  // Toggle for using damage-dependent properties
  params.addParam<bool>("use_damaged_properties", true,
      "Use damage-dependent biot coefficient, porosity, and biot modulus from material properties");

  // Fallback constant values
  params.addParam<Real>("porosity_value", 0.2, "Constant porosity (used if use_damaged_properties=false)");
  params.addParam<Real>("biot_coefficient_value", 0.75, "Constant Biot coefficient (used if use_damaged_properties=false)");

  return params;
}

ElkADPoroDynamicPhaseFieldMaterials::ElkADPoroDynamicPhaseFieldMaterials(const InputParameters & parameters)
  : ADMaterial(parameters),
    _rhos_val(getParam<Real>("rhos_value")),
    _rhof_val(getParam<Real>("rhof_value")),
    _tortosity_val(getParam<Real>("tortosity_value")),
    _viscosity_val(getParam<Real>("viscosity_value")),
    _grain_bulk_modulus(getParam<Real>("grain_bulk_modulus")),
    _fluid_bulk_modulus(getParam<Real>("fluid_bulk_modulus")),
    _use_damaged_properties(getParam<bool>("use_damaged_properties")),
    _porosity_const(getParam<Real>("porosity_value")),
    _biot_coefficient_const(getParam<Real>("biot_coefficient_value")),
    _biot_coefficient_damaged(nullptr),
    _porosity_damaged(nullptr),
    _biot_modulus_damaged(nullptr),
    _effective_perm(getADMaterialProperty<RankTwoTensor>("effective_perm")),
    _elastic_strain(getADMaterialProperty<RankTwoTensor>("elastic_strain")),
    _rhof(declareADProperty<Real>("rhof")),
    _mixture_density(declareADProperty<Real>("mixture_density")),
    _density(declareADProperty<Real>("density")),
    _porosity(declareADProperty<Real>("porosity")),
    _tortosity(declareADProperty<Real>("tortosity")),
    _viscosity(declareADProperty<Real>("viscosity")),
    _biot_modulus(declareADProperty<Real>("biot_modulus")),
    _biot_coefficient(declareADProperty<Real>("biot_coefficient")),
    _permeability(declareADProperty<RankTwoTensor>("permeability")),
    _vol_strain(declareADProperty<Real>("vol_strain"))
{
  if (_use_damaged_properties)
  {
    _biot_coefficient_damaged = &getADMaterialProperty<Real>("biot_coefficient_damaged");
    _porosity_damaged = &getADMaterialProperty<Real>("PorousFlow_porosity_qp_damaged");
    _biot_modulus_damaged = &getADMaterialProperty<Real>("PorousFlow_constant_biot_modulus_qp");
  }
}

void
ElkADPoroDynamicPhaseFieldMaterials::initQpStatefulProperties()
{
  // Fluid density
  _rhof[_qp] = _rhof_val;

  // Get porosity
  ADReal phi = _use_damaged_properties ? (*_porosity_damaged)[_qp] : ADReal(_porosity_const);

  // Mixture density (for reference/output only)
  _mixture_density[_qp] = _rhos_val * (1.0 - phi) + _rhof_val * phi;

  // Constant solid density for ADInertialForce (matches two-field approach for energy consistency)
  _density[_qp] = _rhos_val;

  // Porosity
  _porosity[_qp] = phi;

  // Tortosity
  _tortosity[_qp] = _tortosity_val;

  // Viscosity
  _viscosity[_qp] = _viscosity_val;

  // Biot coefficient
  ADReal alpha = _use_damaged_properties ? (*_biot_coefficient_damaged)[_qp]
                                         : ADReal(_biot_coefficient_const);
  _biot_coefficient[_qp] = alpha;

  // Biot modulus
  if (_use_damaged_properties)
  {
    _biot_modulus[_qp] = (*_biot_modulus_damaged)[_qp];
  }
  else
  {
    // Compute Biot modulus: 1/M = phi/Kf + (alpha - phi)/Ks
    ADReal denom = phi / _fluid_bulk_modulus + (alpha - phi) / _grain_bulk_modulus;
    _biot_modulus[_qp] = 1.0 / std::max(denom, ADReal(1e-20));
  }

  // Permeability (from phase-field elasticity model)
  _permeability[_qp] = _effective_perm[_qp];

  // Volumetric strain
  _vol_strain[_qp] = 0.0;
}

void
ElkADPoroDynamicPhaseFieldMaterials::computeQpProperties()
{
  // Fluid density
  _rhof[_qp] = _rhof_val;

  // Get porosity
  ADReal phi = _use_damaged_properties ? (*_porosity_damaged)[_qp] : ADReal(_porosity_const);

  // Mixture density (for reference/output only, damage-dependent through porosity)
  _mixture_density[_qp] = _rhos_val * (1.0 - phi) + _rhof_val * phi;

  // Constant solid density for ADInertialForce (matches two-field approach for energy consistency)
  _density[_qp] = _rhos_val;

  // Porosity
  _porosity[_qp] = phi;

  // Tortosity
  _tortosity[_qp] = _tortosity_val;

  // Viscosity
  _viscosity[_qp] = _viscosity_val;

  // Biot coefficient
  ADReal alpha = _use_damaged_properties ? (*_biot_coefficient_damaged)[_qp]
                                         : ADReal(_biot_coefficient_const);
  _biot_coefficient[_qp] = alpha;

  // Biot modulus
  if (_use_damaged_properties)
  {
    _biot_modulus[_qp] = (*_biot_modulus_damaged)[_qp];
  }
  else
  {
    // Compute Biot modulus: 1/M = phi/Kf + (alpha - phi)/Ks
    ADReal denom = phi / _fluid_bulk_modulus + (alpha - phi) / _grain_bulk_modulus;
    _biot_modulus[_qp] = 1.0 / std::max(denom, ADReal(1e-20));
  }

  // Permeability (from phase-field elasticity model - already damage-dependent)
  _permeability[_qp] = _effective_perm[_qp];

  // Compute volumetric strain from elastic strain tensor
  _vol_strain[_qp] = _elastic_strain[_qp].trace();
}
