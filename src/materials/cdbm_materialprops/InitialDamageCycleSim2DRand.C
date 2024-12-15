//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "InitialDamageCycleSim2DRand.h"
#include <random>

/**
 *  Created by Chunhui Zhao, Oct 20th, 2024
 *  Material used in defining initial damage profile with exponential decay along the normal direction
 *  This is for 2D only
 */
registerMooseObject("farmsApp", InitialDamageCycleSim2DRand);

InputParameters
InitialDamageCycleSim2DRand::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Material used in defining initial damage profile");
  params.addRequiredParam<Real>("len_of_fault", "length of the fault");
  params.addRequiredParam<Real>("sigma", "decay rate");
  params.addRequiredParam<Real>("peak_val", "peak value of the initial damage");
  params.addParam<bool>("use_background_randalpha", false, "Whether to use random alpha background");
  params.addParam<bool>("use_decay_region", true, "Whether to use decay region");
  params.addCoupledVar("randalpha", "Random alpha auxiliary variable");
  return params;
}

InitialDamageCycleSim2DRand::InitialDamageCycleSim2DRand(const InputParameters & parameters)
  : Material(parameters),
  _initial_damage(declareProperty<Real>("initial_damage")),
  _len_of_fault(getParam<Real>("len_of_fault")),
  _sigma(getParam<Real>("sigma")),
  _peak_val(getParam<Real>("peak_val")),
  _use_background_randalpha(getParam<bool>("use_background_randalpha")),
  _use_decay_region(getParam<bool>("use_decay_region")),
  _randalpha(_use_background_randalpha ? &coupledValue("randalpha") : nullptr)
{
}

void
InitialDamageCycleSim2DRand::initQpStatefulProperties()
{
}

void
InitialDamageCycleSim2DRand::computeQpProperties()
{

  Real x_coord = _q_point[_qp](0); //along the strike direction
  Real y_coord = _q_point[_qp](1); //along the normal direction

  Real alpha_o = 0; // Initialize the damage value
  Real sigma = _sigma; // Decay rate parameter

  if (_use_decay_region){
    // Inside the core region
    if (x_coord > -3000 && x_coord < 3000 && y_coord > -1500 && y_coord < 1500)
    {
      // High damage region
      if (x_coord > -1000 && x_coord < 1000 && y_coord > -200 && y_coord < 200)
      {
        alpha_o = 0.8;
      }
      else
      {
        // Background random alpha
        alpha_o = (*_randalpha)[_qp];
      }
    }
    else
    {
      // Decay region
      Real r = 0.0;

      if (x_coord <= -3000 && y_coord >= -1500 && y_coord <= 1500)
      {
        // Left edge
        r = std::abs(x_coord + 3000);
      }
      else if (x_coord >= 3000 && y_coord >= -1500 && y_coord <= 1500)
      {
        // Right edge
        r = std::abs(x_coord - 3000);
      }
      else if (y_coord <= -1500 && x_coord >= -3000 && x_coord <= 3000)
      {
        // Bottom edge
        r = std::abs(y_coord + 1500);
      }
      else if (y_coord >= 1500 && x_coord >= -3000 && x_coord <= 3000)
      {
        // Top edge
        r = std::abs(y_coord - 1500);
      }
      else if (x_coord >= 3000 && y_coord >= 1500)
      {
        // Top right corner
        r = std::sqrt(std::pow(x_coord - 3000, 2) + std::pow(y_coord - 1500, 2));
      }
      else if (x_coord >= 3000 && y_coord <= -1500)
      {
        // Bottom right corner
        r = std::sqrt(std::pow(x_coord - 3000, 2) + std::pow(y_coord + 1500, 2));
      }
      else if (x_coord <= -3000 && y_coord >= 1500)
      {
        // Top left corner
        r = std::sqrt(std::pow(x_coord + 3000, 2) + std::pow(y_coord - 1500, 2));
      }
      else if (x_coord <= -3000 && y_coord <= -1500)
      {
        // Bottom left corner
        r = std::sqrt(std::pow(x_coord + 3000, 2) + std::pow(y_coord + 1500, 2));
      }

      // Decay from 0.4 to 0 based on distance from the edge
      alpha_o = std::max(0.4 * std::exp(-std::pow(r, 2) / (sigma * sigma)), 0.0);
    }
  }
  else
  {
    // Inside the core region
    if (x_coord > -3000 && x_coord < 3000 && y_coord > -1500 && y_coord < 1500)
    {
      // High damage region
      if (x_coord > -1000 && x_coord < 1000 && y_coord > -200 && y_coord < 200)
      {
        alpha_o = 0.8;
      }
      else
      {
        // Background random alpha
        alpha_o = (*_randalpha)[_qp];
      }
    }
    else
    {
      // Background random alpha
      alpha_o = (*_randalpha)[_qp];
    }

  }

  _initial_damage[_qp] = alpha_o; 
  
}