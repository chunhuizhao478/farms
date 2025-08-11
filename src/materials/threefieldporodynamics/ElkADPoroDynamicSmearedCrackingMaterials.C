//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ElkADPoroDynamicSmearedCrackingMaterials.h"

/**
 *  Created by Chunhui Zhao, Jul 15th, 2024
 *  ADMaterial used in three field poro dynamics smeared cracking simulations
 */
registerMooseObject("farmsApp", ElkADPoroDynamicSmearedCrackingMaterials);

InputParameters
ElkADPoroDynamicSmearedCrackingMaterials::validParams()
{
  InputParameters params = ADMaterial::validParams();
  params.addClassDescription("ADMaterial used in three field poro dynamics simulations");
  params.addRequiredParam<Real>("rhos_value", "solid density value");
  params.addRequiredParam<Real>("rhof_value", "fluid density value");
  params.addRequiredParam<Real>("porosity_value", "porosity value");
  params.addRequiredParam<Real>("tortosity_value", "tortosity value");
  params.addRequiredParam<Real>("viscosity_value", "viscosity value");
  params.addRequiredParam<Real>("bulk_modulus_solid_value", "biot modulus value for solid");
  params.addRequiredParam<Real>("biot_coefficient_value", "biot coefficient value");
  params.addRequiredParam<Real>("bulk_modulus_fluid_value", "biot modulus value for fluid");
  params.addRequiredParam<RealTensorValue>("permeability_value", "permeability value");
  return params;
}

ElkADPoroDynamicSmearedCrackingMaterials::ElkADPoroDynamicSmearedCrackingMaterials(const InputParameters & parameters)
  : ADMaterial(parameters),
    _rhos_val(getParam<Real>("rhos_value")),
    _rhof_val(getParam<Real>("rhof_value")),
    _porosity_val(getParam<Real>("porosity_value")),
    _tortosity_val(getParam<Real>("tortosity_value")),
    _viscosity_val(getParam<Real>("viscosity_value")),
    _bulk_modulus_solid_val(getParam<Real>("bulk_modulus_solid_value")),
    _biot_coefficient_val(getParam<Real>("biot_coefficient_value")),
    _bulk_modulus_fluid_val(getParam<Real>("bulk_modulus_fluid_value")),
    _permeability_val(getParam<RealTensorValue>("permeability_value")),
    _rhof(declareADProperty<Real>("rhof")),
    _rho(declareADProperty<Real>("density")),
    _porosity(declareADProperty<Real>("porosity")),
    _tortosity(declareADProperty<Real>("tortosity")),
    _viscosity(declareADProperty<Real>("viscosity")),
    _biot_modulus(declareADProperty<Real>("biot_modulus")),
    _biot_coefficient(declareADProperty<Real>("biot_coefficient")),
    _permeability(declareADProperty<RealTensorValue>("permeability")),
    _vol_strain(declareADProperty<Real>("vol_strain")),
    _effective_perm(getADMaterialProperty<RealTensorValue>("effective_perm")),
    // _solid_bulk_compliance_damaged(getADMaterialProperty<Real>("solid_bulk_compliance_damaged")),
    _elastic_strain(getADMaterialProperty<RankTwoTensor>("elastic_strain"))
{
}

void
ElkADPoroDynamicSmearedCrackingMaterials::initQpStatefulProperties()
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
  //ADReal biot_coefficient = 1 - (1.0 * _bulk_modulus_solid_skeleton_val) / (1.0 * _bulk_modulus_solid_val);
  ADReal biot_coefficient = _biot_coefficient_val;

  //biot modulus
  //Form adopted in (Mella, 2023)
  //_biot_modulus[_qp] = 1.0 / ( ( _porosity_val / _bulk_modulus_fluid_val ) + ( (biot_coefficient - _porosity_val)/( _bulk_modulus_solid_val ) ) );
  //biot modulus
  // Correct: 1/M = phi/Kf + (b - phi)/Ks
  _biot_modulus[_qp] = 1.0 / ( ( _porosity_val / _bulk_modulus_fluid_val ) +
                               ( (biot_coefficient - _porosity_val) / _bulk_modulus_solid_val ) );

  //biot coefficient
  _biot_coefficient[_qp] = biot_coefficient;

  //permeability
  _permeability[_qp] = _permeability_val;

  //volumetric strain
  _vol_strain[_qp] = 0.0;
}

void
ElkADPoroDynamicSmearedCrackingMaterials::computeQpProperties()
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
  //ADReal biot_coefficient = 1 - (1.0 * _bulk_modulus_solid_skeleton_val) / (1.0 * _bulk_modulus_solid_val);
  ADReal biot_coefficient = _biot_coefficient_val;

  //biot modulus //modify from ElkADComputeSmearedCrackingStress
  //Form adopted in (Mella, 2023)
  //_biot_modulus[_qp] = 1.0 / ( ( _porosity_val / _bulk_modulus_fluid_val ) + ( (biot_coefficient - _porosity_val) * _solid_bulk_compliance_damaged[_qp] ) );
  //Form adopted in PorousFlowConstantBiotModulus in MOOSE
  //_biot_modulus[_qp] = 1.0 / ( ( _porosity_val / _bulk_modulus_fluid_val ) + ( (1 - biot_coefficient)*(biot_coefficient - _porosity_val)/_bulk_modulus_solid_val ) );

  //biot modulus
  // Correct: 1/M = phi/Kf + (b - phi)/Ks
  _biot_modulus[_qp] = 1.0 / ( ( _porosity_val / _bulk_modulus_fluid_val ) +
                               ( (biot_coefficient - _porosity_val) / _bulk_modulus_solid_val ) );

  //biot coefficient
  _biot_coefficient[_qp] = biot_coefficient;

  //permeability //modify from ElkADComputeSmearedCrackingStress
  _permeability[_qp] = _effective_perm[_qp];

  //compute volumetric strain
  _vol_strain[_qp] = _elastic_strain[_qp].trace();

}