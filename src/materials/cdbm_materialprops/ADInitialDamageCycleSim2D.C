//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADInitialDamageCycleSim2D.h"
#include <random>

/**
 *  Created by Chunhui Zhao, Oct 20th, 2024
 *  Material used in defining initial damage profile with exponential decay along the normal direction
 *  This is for 2D only
 */
registerMooseObject("farmsApp", ADInitialDamageCycleSim2D);

InputParameters
ADInitialDamageCycleSim2D::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Material used in defining initial damage profile");
  params.addRequiredParam<std::vector<Real>>("geoparams", "geometry parameters coordinates (fault_ymin, fault_ymax, DB_ymin, DB_ymax)");
  return params;
}

ADInitialDamageCycleSim2D::ADInitialDamageCycleSim2D(const InputParameters & parameters)
  : Material(parameters),
  _initial_damage(declareADProperty<Real>("initial_damage")),
  _geoparams(getParam<std::vector<Real>>("geoparams"))
{
}

void
ADInitialDamageCycleSim2D::initQpStatefulProperties()
{
}

void
ADInitialDamageCycleSim2D::computeQpProperties()
{

  ADReal x_coord = _q_point[_qp](0); //along the strike direction
  ADReal y_coord = _q_point[_qp](1); //along the normal direction

  ADReal alpha_o = 0;

  std::random_device rd;
  std::mt19937 gen(rd());
  std::weibull_distribution<double> wb_distribution(2.0,0.05);

  if (y_coord >= _geoparams[0] and y_coord <= _geoparams[1]){
    alpha_o = 0.7;
  }
  else if (y_coord >= _geoparams[2] and y_coord <= _geoparams[3]){
    alpha_o = 0.0;
  }
  else{
    alpha_o = 0.0;
  }

  _initial_damage[_qp] = alpha_o; 
  
}