//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "IntactPorousSolidProperties.h"

/**
 *  Material used in damage-breakage large deformation formulation, consider full damage evolution equation with diffusion
 *  Created by Chunhui Zhao, Dec 24th, 2024
 */
registerMooseObject("farmsApp", IntactPorousSolidProperties);

InputParameters
IntactPorousSolidProperties::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Material used in three field poro dynamics simulations");
  params.addParam<Real>(                   "lambda_o", "initial lambda constant value");
  params.addParam<Real>(            "shear_modulus_o", "initial shear modulus value");
  params.addParam<Real>(         "permeability_solid_o", "permeability of solid meterial");
  params.addParam<Real>(         "initial_viscosity_fluid", "fluid viscosity");
  params.addParam<Real>(         "solid_bulk_modulus_s", "solid bulk modulus of solid grains");
  params.addParam<Real>(         "fluid_bulk_modulus", "fluid bulk modulus"); 
  params.addParam<Real>(         "porosity_solid_o", "initial prosoity of solid phase"); 
  //input parameters
  return params;
}

IntactPorousSolidProperties::IntactPorousSolidProperties(const InputParameters & parameters)
  : Material(parameters),
    _lambda_o_value(getParam<Real>("lambda_o")),
    _shear_modulus_o_value(getParam<Real>("shear_modulus_o")),
    _permeability_solid_o(getParam<Real>("permeability_solid_o")),
    _solid_bulk_modulus_s(getParam<Real>("solid_bulk_modulus_s")),
    _fluid_bulk_modulus(getParam<Real>("fluid_bulk_modulus")),
    _porosity_solid_o(getParam<Real>("porosity_solid_o")),
    _initial_viscosity_fluid(getParam<Real>("initial_viscosity_fluid")),
    _biot_coeff_eff(declareProperty<Real>("biot_coefficient_effective")),
    _Biot_modulus_eff(declareProperty<Real>("Biot_modulus_effective")),
    _fluid_solid_coupling(declareProperty<Real>("fluid_solid_coupling")),
    _perm_s(declareProperty<Real>("permeability_solid"))
    
{
}

//Rules:See https://github.com/idaholab/moose/discussions/19450
//Only the object that declares the material property can assign values to it.
//Objects can request material properties, gaining read-only access to their values.
//When any object (including the object that declares it) requests the old value of a material property, that property becomes "stateful".
//All stateful material properties must be initialized within the initQpStatefulProperties call. 
//
void
IntactPorousSolidProperties::computeQpProperties()
{

  _biot_coeff_eff[_qp] = 1 - (_lambda_o_value +  0.666667 * _shear_modulus_o_value) / _solid_bulk_modulus_s;

  _fluid_solid_coupling[_qp] = _biot_coeff_eff[_qp];

  Real denemnator = (_biot_coeff_eff[_qp] - _porosity_solid_o) * _solid_bulk_modulus_s + _porosity_solid_o * _fluid_bulk_modulus;

  _Biot_modulus_eff[_qp] = _fluid_bulk_modulus * _solid_bulk_modulus_s / denemnator;

  _perm_s[_qp] = _permeability_solid_o / _initial_viscosity_fluid;
  
}

