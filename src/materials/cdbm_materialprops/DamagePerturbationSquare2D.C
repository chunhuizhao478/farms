//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DamagePerturbationSquare2D.h"

/**
 *  Created by Chunhui Zhao, Aug 27th, 2024
 *  Material used in Create Time Dependent Damage Perturbation in the Dynamic Solve
 */
registerMooseObject("farmsApp", DamagePerturbationSquare2D);

InputParameters
DamagePerturbationSquare2D::validParams()
{
  InputParameters params = ADMaterial::validParams();
  params.addClassDescription("Material used in Create Time Dependent Damage Perturbation in the Dynamic Solve");
  params.addRequiredParam<std::vector<Real>>("nucl_center", "nucleation center (x,y,z)");
  params.addRequiredParam<Real>("e_damage","the peak damage value apply region normal to the fault (exponential decay)");
  params.addRequiredParam<Real>("thickness","the standard deviation used in apply damage value normal to the fault (exponential decay)");
  params.addRequiredParam<Real>("length","the standard deviation used in apply damage value normal to the fault (exponential decay)");
  params.addRequiredParam<Real>("duration","duration to reach peak damage");
  params.addRequiredParam<Real>("sigma","the standard deviation used in apply damage value normal to the fault (exponential decay)");
  return params;
}

DamagePerturbationSquare2D::DamagePerturbationSquare2D(const InputParameters & parameters)
  : Material(parameters),
  _damage_perturbation(declareProperty<Real>("damage_perturb")),
  _damage_perturbation_old(getMaterialPropertyOldByName<Real>("damage_perturb")),
  _nucl_center(getParam<std::vector<Real>>("nucl_center")),
  _peak_damage(getParam<Real>("e_damage")),
  _thickness(getParam<Real>("thickness")),
  _length(getParam<Real>("length")),
  _duration(getParam<Real>("duration")),
  _sigma(getParam<Real>("sigma"))
{
}

void
DamagePerturbationSquare2D::initQpStatefulProperties()
{
  _damage_perturbation[_qp] = 0.0;
}

void
DamagePerturbationSquare2D::computeQpProperties()
{

  //Get coordinates
  //here no rotation is applied yet
  Real xcoord = _q_point[_qp](0); //strike
  Real ycoord = _q_point[_qp](1); //normal

  //Get damage increments
  Real damage_inc = _peak_damage / (_duration / _dt);

  //Assign initial damage perturbation
  Real dalpha = 0.0;
  if ( _t <= _duration ){
    if ( (xcoord >= _nucl_center[0] - _length / 2.0) && (xcoord <= _nucl_center[0] + _length / 2.0) && (ycoord >= _nucl_center[1] - _thickness / 2.0) && (ycoord <= _nucl_center[1] + _thickness / 2.0) ){
      dalpha = _damage_perturbation_old[_qp] + damage_inc*std::exp(-(xcoord*xcoord)/(2.0*_sigma*_sigma));
      //dalpha = _peak_damage;
    }
    else{
      dalpha = _damage_perturbation_old[_qp];
    }
  }
  else{
    dalpha = _damage_perturbation_old[_qp];
  }
  _damage_perturbation[_qp] = dalpha;
}