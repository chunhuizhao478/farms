//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "InitialDamageCycleSim3D.h"
#include <random>

/**
 *  Created by Chunhui Zhao, Nov 28th, 2024
 *  Material used in defining initial damage profile with exponential decay along the normal direction
 *  This is for 3D only
 */
registerMooseObject("farmsApp", InitialDamageCycleSim3D);

InputParameters
InitialDamageCycleSim3D::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Material used in defining initial damage profile");
  params.addRequiredParam<Real>("len_of_fault", "length of the fault");
  params.addRequiredParam<Real>("len_along_dip", "length along fault dip direction");
  params.addRequiredParam<Real>("sigma", "decay rate");
  params.addRequiredParam<Real>("peak_val", "peak value of the initial damage");
  params.addRequiredParam<std::vector<Real>>("nucl_center", "nucleation center (x,y,z)");
  return params;
}

InitialDamageCycleSim3D::InitialDamageCycleSim3D(const InputParameters & parameters)
  : Material(parameters),
  _initial_damage(declareProperty<Real>("initial_damage")),
  _len_of_fault(getParam<Real>("len_of_fault")),
  _len_along_dip(getParam<Real>("len_along_dip")),
  _sigma(getParam<Real>("sigma")),
  _peak_val(getParam<Real>("peak_val")),
  _nucl_center(getParam<std::vector<Real>>("nucl_center"))
{
}

void
InitialDamageCycleSim3D::initQpStatefulProperties()
{
}

void
InitialDamageCycleSim3D::computeQpProperties()
{

  Real x_coord = _q_point[_qp](0); // along strike
  Real y_coord = _q_point[_qp](1); // along normal
  Real z_coord = _q_point[_qp](2); // along dip

  Real alpha_o = 0.0;
  Real r = 0.0;
  Real sigma = _sigma;

  // Check if point is within fault plane
  if (x_coord > -0.5*_len_of_fault && 
      x_coord < 0.5*_len_of_fault && 
      z_coord < 0 && 
      z_coord > -_len_along_dip)
  {
    // Distance from fault plane
    r = std::abs(y_coord);
  }
  else
  {
    // Distance to nearest fault boundary
    Real dist_x = std::min(std::abs(x_coord - 0.5*_len_of_fault), 
                          std::abs(x_coord + 0.5*_len_of_fault));
    Real dist_y = std::abs(y_coord);
    Real dist_z = std::min(std::abs(z_coord), 
                          std::abs(z_coord + _len_along_dip));
    
    // Compute minimum distance to boundaries
    r = std::sqrt(dist_x*dist_x + dist_y*dist_y + dist_z*dist_z);
  }

  // Apply Gaussian decay
  alpha_o = _peak_val * std::exp(-1.0 * (r*r)/(sigma*sigma));
  
  // Ensure non-negative values
  _initial_damage[_qp] = std::max(alpha_o, 0.0);
  
}