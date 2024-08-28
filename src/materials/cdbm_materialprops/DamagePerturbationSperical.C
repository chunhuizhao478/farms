//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DamagePerturbationSperical.h"

/**
 *  Created by Chunhui Zhao, Aug 27th, 2024
 *  Material used in Create Time Dependent Damage Perturbation in the Dynamic Solve
 */
registerMooseObject("farmsApp", DamagePerturbationSperical);

InputParameters
DamagePerturbationSperical::validParams()
{
  InputParameters params = ADMaterial::validParams();
  params.addClassDescription("Material used in Create Time Dependent Damage Perturbation in the Dynamic Solve");
  params.addRequiredParam<std::vector<Real>>("nucl_center", "nucleation center (x,y,z)");
  params.addRequiredParam<Real>("e_damage","the peak damage value apply region normal to the fault (exponential decay)");
  params.addRequiredParam<Real>("e_sigma","the standard deviation used in apply damage value normal to the fault (exponential decay)");
  params.addRequiredParam<Real>("duration","duration to reach peak damage");
  return params;
}

DamagePerturbationSperical::DamagePerturbationSperical(const InputParameters & parameters)
  : Material(parameters),
  _damage_perturbation(declareProperty<Real>("damage_perturbation")),
  _damage_perturbation_old(getMaterialPropertyOldByName<Real>("damage_perturbation")),
  _initial_damage(getMaterialProperty<Real>("initial_damage")),
  _nucl_center(getParam<std::vector<Real>>("nucl_center")),
  _peak_damage(getParam<Real>("e_damage")),
  _sigma(getParam<Real>("e_sigma")),
  _duration(getParam<Real>("duration"))
{
  //in case I'm stupid
  if (_nucl_center.size() != 3){
    mooseError("fault plane parameter must accepts 3 numbers!");
  }
}

void
DamagePerturbationSperical::initQpStatefulProperties()
{
  _damage_perturbation[_qp] = 0.0;
}

void
DamagePerturbationSperical::computeQpProperties()
{

  //Get coordinates
  //here no rotation is applied yet
  Real xcoord = _q_point[_qp](0); //strike
  Real ycoord = _q_point[_qp](1); //dip
  Real zcoord = _q_point[_qp](2); //normal

  //Get damage increments
  Real damage_inc = _peak_damage / (_duration / _dt);

  //Assign initial damage perturbation
  Real dalpha = 0.0;
  if ( _t <= _duration ){
    Real r = std::sqrt(pow(xcoord - _nucl_center[0],2) + pow(ycoord - _nucl_center[1],2) + pow(zcoord - _nucl_center[2],2));
    dalpha = _damage_perturbation_old[_qp] + std::max(damage_inc * exp(-1.0*(std::pow(r,2))/(_sigma*_sigma)),0.0);
  }
  else{
    dalpha = _damage_perturbation_old[_qp];
  }
  _damage_perturbation[_qp] = dalpha;
}