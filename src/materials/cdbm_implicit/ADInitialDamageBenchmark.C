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
  if (_fault_plane.size() != 4){ 
    mooseError("fault plane parameter must accept 4 numbers!");
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
  
  //restrict ourselves to fault plane only, but extend damage along normal direction with exponential decaying
  if ((ycoord >= _fault_plane[2]) && (ycoord <= _fault_plane[3]))
  { 
    //set a constant damage value for nucleation region
    if ( (xcoord >= _nucl_center[0] - 0.5 * _nucl_distance) && (xcoord <= _nucl_center[0] + 0.5 * _nucl_distance) && (ycoord >= _nucl_center[1] - 0.5 * _nucl_distance) && (ycoord <= _nucl_center[1] + 0.5 * _nucl_distance) && (zcoord >= _nucl_center[2] - 0.5 * _nucl_thickness) && (zcoord <= _nucl_center[2] + 0.5 * _nucl_thickness) ){
      alpha_o = _nucl_damage;
    } 
    else{ //set a exponential decaying around the fault
      //end of fault left
      if (xcoord < _fault_plane[0]){
        Real r = std::sqrt(pow(xcoord - _fault_plane[0],2) + pow(zcoord - 0,2));
        alpha_o = std::max(_peak_damage * exp(-1.0*(std::pow(r,2))/(_sigma*_sigma)),0.0);
      }
      //end of fault right
      else if (xcoord > _fault_plane[1]){
        Real r = std::sqrt(pow(xcoord - _fault_plane[1],2) + pow(zcoord - 0,2));
        alpha_o = std::max(_peak_damage * exp(-1.0*(std::pow(r,2))/(_sigma*_sigma)),0.0);
      }
      //along the fault
      else{
        alpha_o = std::max(_peak_damage * exp(-1.0*(zcoord*zcoord)/(_sigma*_sigma)),0.0);
      }
    }
  }
  _initial_damage[_qp] = alpha_o; 
}