//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ElkADPoroDynamicMaterials.h"

/**
 *  Created by Chunhui Zhao, Jul 6th, 2024
 *  ADMaterial used in three field poro dynamics simulations
 */
registerMooseObject("farmsApp", ElkADPoroDynamicMaterials);

InputParameters
ElkADPoroDynamicMaterials::validParams()
{
  InputParameters params = ADMaterial::validParams();
  params.addClassDescription("ADMaterial used in three field poro dynamics simulations");
  params.addRequiredParam<Real>("rhos_value", "solid density value");
  params.addRequiredParam<Real>("rhof_value", "fluid density value");
  params.addRequiredParam<Real>("porosity_value", "porosity value");
  params.addRequiredParam<Real>("tortosity_value", "tortosity value");
  params.addRequiredParam<Real>("viscosity_value", "viscosity value");
  params.addRequiredParam<Real>("bulk_modulus_solid_skeleton_value", "biot modulus value for solid skeleton");
  params.addRequiredParam<Real>("bulk_modulus_solid_value", "biot modulus value for solid");
  params.addRequiredParam<Real>("bulk_modulus_fluid_value", "biot modulus value for fluid");
  params.addRequiredParam<Real>("permeability_value", "permeability value");
  return params;
}

ElkADPoroDynamicMaterials::ElkADPoroDynamicMaterials(const InputParameters & parameters)
  : ADMaterial(parameters),
    _rhos_val(getParam<Real>("rhos_value")),
    _rhof_val(getParam<Real>("rhof_value")),
    _porosity_val(getParam<Real>("porosity_value")),
    _tortosity_val(getParam<Real>("tortosity_value")),
    _viscosity_val(getParam<Real>("viscosity_value")),
    _bulk_modulus_solid_skeleton_val(getParam<Real>("bulk_modulus_solid_skeleton_value")),
    _bulk_modulus_solid_val(getParam<Real>("bulk_modulus_solid_value")),
    _bulk_modulus_fluid_val(getParam<Real>("bulk_modulus_fluid_value")),
    _permeability_val(getParam<Real>("permeability_value")),
    _rhof(declareADProperty<Real>("rhof")),
    _rho(declareADProperty<Real>("density")),
    _porosity(declareADProperty<Real>("porosity")),
    _tortosity(declareADProperty<Real>("tortosity")),
    _viscosity(declareADProperty<Real>("viscosity")),
    _biot_modulus(declareADProperty<Real>("biot_modulus")),
    _biot_coefficient(declareADProperty<Real>("biot_coefficient")),
    _permeability(declareADProperty<Real>("permeability"))
{
}

void
ElkADPoroDynamicMaterials::computeQpProperties()
{
  //fluid density
  _rhof[_qp] = _rhof_val;

  //density
  _rho[_qp] = _rhos_val * ( 1 - _porosity_val ) + _rhof_val * _porosity_val;
  
  //porosity
  _porosity[_qp] = _porosity_val;
  
  //tortosity
  _tortosity[_qp] = _tortosity_val;
  
  //viscosity
  _viscosity[_qp] = _viscosity_val;
  
  //* biot coefficient
  ADReal biot_coefficient = 1 - (1.0 * _bulk_modulus_solid_skeleton_val) / (1.0 * _bulk_modulus_solid_val);

  //biot modulus
  _biot_modulus[_qp] = 1.0 / ( ( _porosity_val / _bulk_modulus_fluid_val ) + ( (biot_coefficient - _porosity_val)/( _bulk_modulus_solid_val ) ) );
  
  //biot coefficient
  _biot_coefficient[_qp] = biot_coefficient;

  //permeability
  _permeability[_qp] = _permeability_val;
}