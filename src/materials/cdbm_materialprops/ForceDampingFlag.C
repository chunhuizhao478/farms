//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ForceDampingFlag.h"

/**
 *  Created by Chunhui Zhao, Nov 19th, 2024
 *  Material used in get the flag for force damping based on the velocity
 */
registerMooseObject("farmsApp", ForceDampingFlag);

InputParameters
ForceDampingFlag::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Material used in defining flag for force damping based on velocity");
  params.addRequiredParam<Real>("vel_maximum_threshold", "The maximum velocity threshold for force damping");
  params.addRequiredParam<Real>("vel_minimum_threshold", "The minimum velocity threshold for force damping");
  params.addRequiredParam<PostprocessorName>("max_vel_x", "Postprocessor for maximum velocity");  // Add postprocessor param
  params.addRequiredParam<PostprocessorName>("max_vel_y", "Postprocessor for maximum velocity");
  return params;
}

ForceDampingFlag::ForceDampingFlag(const InputParameters & parameters)
  : Material(parameters),
  _flag(declareProperty<Real>("flag")),
  _flag_old(getMaterialProperty<Real>("flag")),
  _vel_maximum_threshold(getParam<Real>("vel_maximum_threshold")),
  _vel_minimum_threshold(getParam<Real>("vel_minimum_threshold")),
  _max_vel_x(getPostprocessorValue("max_vel_x")),  // Get postprocessor value
  _max_vel_y(getPostprocessorValue("max_vel_y"))
{
}

void
ForceDampingFlag::initQpStatefulProperties()
{
  _flag[_qp] = 0;
}

void
ForceDampingFlag::computeQpProperties()
{
  Real max_vel_mag = sqrt(_max_vel_x * _max_vel_x + _max_vel_y * _max_vel_y);
  // 
  if ( max_vel_mag > _vel_maximum_threshold )
  {
    _flag[_qp] = 1;
  }
  else if ( max_vel_mag < _vel_minimum_threshold )
  {
    _flag[_qp] = 0;
  }
  else
  {
    _flag[_qp] = _flag_old[_qp];
  }

}