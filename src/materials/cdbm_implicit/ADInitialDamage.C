//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADInitialDamage.h"

/**
 *  Created by Chunhui Zhao, Jul 18th, 2024
 *  ADMaterial used in defining initial damage profile
 */
registerMooseObject("farmsApp", ADInitialDamage);

InputParameters
ADInitialDamage::validParams()
{
  InputParameters params = ADMaterial::validParams();
  params.addClassDescription("ADMaterial used in defining initial damage profile");
  return params;
}

ADInitialDamage::ADInitialDamage(const InputParameters & parameters)
  : ADMaterial(parameters),
  _initial_damage(declareADProperty<Real>("initial_damage")),
  _initial_damage_old(getMaterialPropertyOldByName<Real>("initial_damage"))
{
}

void
ADInitialDamage::initQpStatefulProperties()
{

  Real xcoord = _q_point[_qp](0);
  Real ycoord = _q_point[_qp](1);
  Real zcoord = _q_point[_qp](2);

  //
  //note: hardcode coordinates
  Real alpha_o = 0.0;
  if (xcoord >= -30 && xcoord <= 30 && ycoord >= -5 && ycoord <= 5 && zcoord >= -0.5 && zcoord <= 0.5){
    if ( xcoord >= -25 && xcoord <= -20 && ycoord >= -2.5 && ycoord <= 2.5 ){
      alpha_o = 0.8;
    }
    else{
      alpha_o = 0.7;
    }
  }
  else{
    alpha_o = 0.0;
  }  

  //
  _initial_damage[_qp] = alpha_o;

}

void
ADInitialDamage::computeQpProperties()
{
  _initial_damage[_qp] = _initial_damage_old[_qp];
}