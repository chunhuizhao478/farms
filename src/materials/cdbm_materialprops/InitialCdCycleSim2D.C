//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "InitialCdCycleSim2D.h"
#include <random>

/**
 *  Created by Chunhui Zhao, Jan 29th, 2025
 *  ADMaterial used in defining initial Cd profile with exponential decay along the normal direction
 *  This is for 2D only
 */
registerMooseObject("farmsApp", InitialCdCycleSim2D);

InputParameters
InitialCdCycleSim2D::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Material used in defining initial damage profile");
  params.addRequiredParam<Real>("len_of_fault", "length of the fault");
  params.addRequiredParam<Real>("sigma", "decay rate");
  params.addRequiredParam<Real>("peak_val", "peak value of the initial damage");
  params.addParam<bool>("use_background_randalpha", false, "Whether to use random alpha background");
  params.addCoupledVar("randalpha", "Random alpha auxiliary variable");
  return params;
}

InitialCdCycleSim2D::InitialCdCycleSim2D(const InputParameters & parameters)
  : Material(parameters),
  _initial_damage(declareProperty<Real>("initial_cd")),
  _len_of_fault(getParam<Real>("len_of_fault")),
  _sigma(getParam<Real>("sigma")),
  _peak_val(getParam<Real>("peak_val")),
  _use_background_randalpha(getParam<bool>("use_background_randalpha")),
  _randalpha(_use_background_randalpha ? &coupledValue("randalpha") : nullptr)
{
}

void
InitialCdCycleSim2D::initQpStatefulProperties()
{
}

void
InitialCdCycleSim2D::computeQpProperties()
{

  Real x_coord = _q_point[_qp](0); //along the strike direction
  Real y_coord = _q_point[_qp](1); //along the normal direction

  Real alpha_o = 0;

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

  if (_use_background_randalpha)
  {
    alpha_o = std::max(alpha_o, (*_randalpha)[_qp]);
  }

  _initial_damage[_qp] = alpha_o; 
  
}