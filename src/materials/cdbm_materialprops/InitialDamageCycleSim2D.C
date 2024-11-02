//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "InitialDamageCycleSim2D.h"
#include <random>

/**
 *  Created by Chunhui Zhao, Oct 20th, 2024
 *  Material used in defining initial damage profile with exponential decay along the normal direction
 *  This is for 2D only
 */
registerMooseObject("farmsApp", InitialDamageCycleSim2D);

InputParameters
InitialDamageCycleSim2D::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Material used in defining initial damage profile");
  return params;
}

InitialDamageCycleSim2D::InitialDamageCycleSim2D(const InputParameters & parameters)
  : Material(parameters),
  _initial_damage(declareProperty<Real>("initial_damage"))
{
}

void
InitialDamageCycleSim2D::initQpStatefulProperties()
{
}

void
InitialDamageCycleSim2D::computeQpProperties()
{

  Real x_coord = _q_point[_qp](0); //along the strike direction
  Real y_coord = _q_point[_qp](1); //along the normal direction

  Real alpha_o = 0;

  // std::random_device rd;
  // std::mt19937 gen(rd());
  // std::weibull_distribution<double> wb_distribution(2.0,0.05);

  // Real value = 1.0;
  // while( value > 0.1 ){
  //   value = wb_distribution(gen);
  // }

  Real r = 0.0;
  Real sigma = 2e2;
  if (x_coord > -10000 and x_coord < 10000){
    r = y_coord;
    alpha_o = std::max(0.7 * std::exp(-1.0*(std::pow(r,2))/(sigma*sigma)),0.0);
  }
  else if (x_coord <= -10000){
    r = std::sqrt((y_coord - 0) * (y_coord - 0) + (x_coord - (-10000)) * (x_coord - (-10000)));
    alpha_o = std::max(0.7 * std::exp(-1.0*(std::pow(r,2))/(sigma*sigma)),0.0);
  }
  else if (x_coord >= 10000){
    r = std::sqrt((y_coord - 0) * (y_coord - 0) + (x_coord - (10000)) * (x_coord - (10000)));
    alpha_o = std::max(0.7 * std::exp(-1.0*(std::pow(r,2))/(sigma*sigma)),0.0);
  }

  _initial_damage[_qp] = alpha_o; 
  
}