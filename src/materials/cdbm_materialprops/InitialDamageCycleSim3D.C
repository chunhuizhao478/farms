//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "InitialDamageCycleSim3D.h"
#include <random>

/**
 *  Created by Chunhui Zhao, Nov 28th, 2024
 *  Material used in defining initial damage profile with exponential decay along the normal direction
 *  This is for 3D only
 */
registerMooseObject("farmsApp", InitialDamageCycleSim3D);

InputParameters
InitialDamageCycleSim3D::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Material used in defining initial damage profile");
  params.addRequiredParam<Real>("len_of_fault", "length of the fault");
  params.addRequiredParam<Real>("len_along_dip", "length along fault dip direction");
  params.addRequiredParam<Real>("sigma", "decay rate");
  params.addRequiredParam<Real>("peak_val", "peak value of the initial damage");
  params.addRequiredParam<std::vector<Real>>("nucl_center", "nucleation center (x,y,z)");
  return params;
}

InitialDamageCycleSim3D::InitialDamageCycleSim3D(const InputParameters & parameters)
  : Material(parameters),
  _initial_damage(declareProperty<Real>("initial_damage")),
  _len_of_fault(getParam<Real>("len_of_fault")),
  _len_along_dip(getParam<Real>("len_along_dip")),
  _sigma(getParam<Real>("sigma")),
  _peak_val(getParam<Real>("peak_val")),
  _nucl_center(getParam<std::vector<Real>>("nucl_center"))
{
}

void
InitialDamageCycleSim3D::initQpStatefulProperties()
{
}

void
InitialDamageCycleSim3D::computeQpProperties()
{

  Real x_coord = _q_point[_qp](0); //along the strike direction
  Real y_coord = _q_point[_qp](1); //along the normal direction
  Real z_coord = _q_point[_qp](2); //along the dip direction

  Real alpha_o = 0;

  Real r = 0.0;
  Real sigma = _sigma;
  if (x_coord > -0.5*_len_of_fault and x_coord < 0.5*_len_of_fault and z_coord < 0 and z_coord > -_len_along_dip){
    if (y_coord >= 0 -0.5 * 400 && y_coord <= 0 + 0.5 * 400 && (x_coord >= 0 - 400 / 2.0) && (x_coord <= 0 + 400 / 2.0) && (z_coord >= -6000 - 400 / 2.0) && (z_coord <= -6000 + 400 / 2.0)){ //set high damage strip
        alpha_o = 0.9;
    }
    else{
      r = y_coord;
      alpha_o = std::max(_peak_val * std::exp(-1.0*(std::pow(r,2))/(sigma*sigma)),0.0);
    }
  }
  _initial_damage[_qp] = alpha_o; 
  
}