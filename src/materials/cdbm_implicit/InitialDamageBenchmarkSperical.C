//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "InitialDamageBenchmarkSperical.h"

/**
 *  Created by Chunhui Zhao, Aug 14th, 2024
 *  ADMaterial used in defining initial damage profile with exponential decay along the normal direction
 */
registerMooseObject("farmsApp", InitialDamageBenchmarkSperical);

InputParameters
InitialDamageBenchmarkSperical::validParams()
{
  InputParameters params = ADMaterial::validParams();
  params.addClassDescription("ADMaterial used in defining initial damage profile");
  params.addRequiredParam<std::vector<Real>>("nucl_center", "nucleation center (x,y,z)");
  params.addRequiredParam<Real>("e_damage","the peak damage value apply region normal to the fault (exponential decay)");
  params.addRequiredParam<Real>("e_sigma","the standard deviation used in apply damage value normal to the fault (exponential decay)");
  return params;
}

InitialDamageBenchmarkSperical::InitialDamageBenchmarkSperical(const InputParameters & parameters)
  : Material(parameters),
  _initial_damage(declareProperty<Real>("initial_damage")),
  _nucl_center(getParam<std::vector<Real>>("nucl_center")),
  _peak_damage(getParam<Real>("e_damage")),
  _sigma(getParam<Real>("e_sigma"))
{
  //in case I'm stupid
  if (_nucl_center.size() != 3){
    mooseError("fault plane parameter must accepts 3 numbers!");
  }
}

void
InitialDamageBenchmarkSperical::initQpStatefulProperties()
{
}

void
InitialDamageBenchmarkSperical::computeQpProperties()
{

  //Get coordinates
  //here no rotation is applied yet
  Real xcoord = _q_point[_qp](0); //strike
  Real ycoord = _q_point[_qp](1); //dip
  Real zcoord = _q_point[_qp](2); //normal

  //Assign initial damage
  Real alpha_o = 0.0;
  Real r = std::sqrt(pow(xcoord - _nucl_center[0],2) + pow(ycoord - _nucl_center[1],2) + pow(zcoord - _nucl_center[2],2));
  alpha_o = std::max(_peak_damage * exp(-1.0*(std::pow(r,2))/(_sigma*_sigma)),0.0);
  _initial_damage[_qp] = alpha_o; 
}