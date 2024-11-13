//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "InitialDamageCycleSim2DDebug.h"
#include <random>

/**
 *  Created by Chunhui Zhao, Oct 20th, 2024
 *  Material used in defining initial damage profile with exponential decay along the normal direction
 *  This is for 2D only
 */
registerMooseObject("farmsApp", InitialDamageCycleSim2DDebug);

InputParameters
InitialDamageCycleSim2DDebug::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Material used in defining initial damage profile");
  params.addRequiredParam<Real>("len_of_fault", "length of the fault");
  params.addRequiredParam<Real>("sigma", "decay rate");
  params.addRequiredParam<Real>("peak_val", "peak value of the initial damage");
  return params;
}

InitialDamageCycleSim2DDebug::InitialDamageCycleSim2DDebug(const InputParameters & parameters)
  : Material(parameters),
  _initial_damage(declareProperty<Real>("initial_damage")),
  _len_of_fault(getParam<Real>("len_of_fault")),
  _sigma(getParam<Real>("sigma")),
  _peak_val(getParam<Real>("peak_val"))
{
}

void
InitialDamageCycleSim2DDebug::initQpStatefulProperties()
{
}

void
InitialDamageCycleSim2DDebug::computeQpProperties()
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
  Real sigma = _sigma;
  if (x_coord > -0.5*_len_of_fault and x_coord < 0.5*_len_of_fault ){
    r = y_coord;
    alpha_o = std::max(_peak_val * std::exp(-1.0*(std::pow(r,2))/(sigma*sigma)),0.0);
  }
  else if (x_coord <= -0.5*_len_of_fault ){
    r = std::sqrt((y_coord - 0) * (y_coord - 0) + (x_coord - (-0.5*_len_of_fault )) * (x_coord - (-0.5*_len_of_fault )));
    alpha_o = std::max(_peak_val * std::exp(-1.0*(std::pow(r,2))/(sigma*sigma)),0.0);
  }
  else if (x_coord >= 0.5*_len_of_fault ){
    r = std::sqrt((y_coord - 0) * (y_coord - 0) + (x_coord - (0.5*_len_of_fault)) * (x_coord - (0.5*_len_of_fault )));
    alpha_o = std::max(_peak_val * std::exp(-1.0*(std::pow(r,2))/(sigma*sigma)),0.0);
  }

  _initial_damage[_qp] = alpha_o; 
  
}