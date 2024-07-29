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
    if ( xcoord >= -26 && xcoord <= -22 && ycoord >= -2 && ycoord <= 2 ){
      Real center_x = -24; Real center_y = 0;
      Real L = 4.0;
      Real A = 1.0; Real sigma = 0.8;
      Real D = std::sqrt((center_x-xcoord)*(center_x-xcoord)+(center_y-ycoord)*(center_y-ycoord));
      Real V = A * std::exp(-D*D/(2*sigma*sigma));
      Real Vmax = A;
      Real d_max = (L * std::sqrt(2))/2;
      Real Vmin = A * std::exp(-(d_max*d_max)/(2*sigma*sigma));
      Real Vnormalized = (V - Vmin) / (Vmax - Vmin);
      alpha_o = Vnormalized * ( 0.8 - 0.7 ) + 0.7;
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

  Real xcoord = _q_point[_qp](0);
  Real ycoord = _q_point[_qp](1);
  Real zcoord = _q_point[_qp](2);

  //
  //note: hardcode coordinates
  Real alpha_o = 0.0;
  if (xcoord >= -30 && xcoord <= 30 && ycoord >= -5 && ycoord <= 5 && zcoord >= -0.5 && zcoord <= 0.5){
    if ( xcoord >= -26 && xcoord <= -22 && ycoord >= -2 && ycoord <= 2 ){
      Real center_x = -24; Real center_y = 0;
      Real L = 4.0;
      Real A = 1.0; Real sigma = 0.8;
      Real D = std::sqrt((center_x-xcoord)*(center_x-xcoord)+(center_y-ycoord)*(center_y-ycoord));
      Real V = A * std::exp(-D*D/(2*sigma*sigma));
      Real Vmax = A;
      Real d_max = (L * std::sqrt(2))/2;
      Real Vmin = A * std::exp(-(d_max*d_max)/(2*sigma*sigma));
      Real Vnormalized = (V - Vmin) / (Vmax - Vmin);
      alpha_o = Vnormalized * ( 0.8 - 0.7 ) + 0.7;
    }
    else{
      alpha_o = 0.7;
    }
  }
  else{
    alpha_o = 0.0;
  }  

  _initial_damage[_qp] = alpha_o;
}