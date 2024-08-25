//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADInitialDamageBenchmark.h"

/**
 *  Created by Chunhui Zhao, Aug 14th, 2024
 *  ADMaterial used in defining initial damage profile with exponential decay along the normal direction
 */
registerMooseObject("farmsApp", ADInitialDamageBenchmark);

InputParameters
ADInitialDamageBenchmark::validParams()
{
  InputParameters params = ADMaterial::validParams();
  params.addClassDescription("ADMaterial used in defining initial damage profile");
  params.addRequiredParam<std::vector<Real>>("nucl_center", "nucleation center (x,y,z)");
  params.addRequiredParam<std::vector<Real>>("fault_plane", "fault plane coordinates (xmin, xmax, ymin, ymax)");
  params.addRequiredParam<Real>("nucl_distance","nucleation distance away from nucleation center along the fault, assuming box shape");
  params.addRequiredParam<Real>("nucl_thickness","nucleation thickness normal to the nucleation center, assume box shape");
  params.addRequiredParam<Real>("nucl_damage","the damage value apply within the nucleation region");
  params.addRequiredParam<Real>("e_damage","the peak damage value apply region normal to the fault (exponential decay)");
  params.addRequiredParam<Real>("e_sigma","the standard deviation used in apply damage value normal to the fault (exponential decay)");
  return params;
}

ADInitialDamageBenchmark::ADInitialDamageBenchmark(const InputParameters & parameters)
  : ADMaterial(parameters),
  _initial_damage(declareADProperty<Real>("initial_damage")),
  _nucl_center(getParam<std::vector<Real>>("nucl_center")),
  _fault_plane(getParam<std::vector<Real>>("fault_plane")),
  _nucl_distance(getParam<Real>("nucl_distance")),
  _nucl_thickness(getParam<Real>("nucl_thickness")),
  _nucl_damage(getParam<Real>("nucl_damage")),
  _peak_damage(getParam<Real>("e_damage")),
  _sigma(getParam<Real>("e_sigma"))
{
  //in case I'm stupid
  if (_fault_plane.size() != 6){ 
    mooseError("fault plane parameter must accept 6 numbers!");
  }
  if (_nucl_center.size() != 3){
    mooseError("fault plane parameter must accepts 3 numbers!");
  }
}

void
ADInitialDamageBenchmark::initQpStatefulProperties()
{
}

void
ADInitialDamageBenchmark::computeQpProperties()
{

  //Get coordinates
  //here no rotation is applied yet
  Real xcoord = _q_point[_qp](0); //strike
  Real ycoord = _q_point[_qp](1); //dip
  Real zcoord = _q_point[_qp](2); //normal

  //Assign initial damage
  Real alpha_o = 0.0;
  
  Real half_width_x = ( _fault_plane[1] - _fault_plane[0] ) / 2.0;
  Real half_width_y = ( _fault_plane[3] - _fault_plane[2] ) / 2.0;
  Real half_width_z = ( _fault_plane[5] - _fault_plane[4] ) / 2.0;

  // Calculate the distance from the closest point on the boundary of the box
  double dx = std::max(std::abs(xcoord - _nucl_center[0]) - half_width_x, 0.0);
  double dy = std::max(std::abs(ycoord - _nucl_center[1]) - half_width_y, 0.0);
  double dz = std::max(std::abs(zcoord - _nucl_center[2]) - half_width_z, 0.0);
  double r = std::sqrt(dx * dx + dy * dy + dz * dz);

  // If inside the box, return 0.6
  if (dx == 0.0 && dy == 0.0 && dz == 0.0) {
    if ( (xcoord >= _nucl_center[0] - 0.5 * _nucl_distance) && (xcoord <= _nucl_center[0] + 0.5 * _nucl_distance) && (ycoord >= _nucl_center[1] - 0.5 * _nucl_distance) && (ycoord <= _nucl_center[1] + 0.5 * _nucl_distance) && (zcoord >= _nucl_center[2] - 0.5 * _nucl_thickness) && (zcoord <= _nucl_center[2] + 0.5 * _nucl_thickness) ){
      alpha_o = _nucl_damage;
    } 
    else{
      alpha_o = _peak_damage;
    }
  }
  else{
    alpha_o = _peak_damage * std::exp(-1.0 * std::pow(r, 2) / (_sigma * _sigma));
  }

  _initial_damage[_qp] = alpha_o; 
}