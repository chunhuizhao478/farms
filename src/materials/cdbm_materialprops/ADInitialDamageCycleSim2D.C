//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADInitialDamageCycleSim2D.h"
#include <random>

/**
 *  Created by Chunhui Zhao, Oct 20th, 2024
 *  Material used in defining initial damage profile with exponential decay along the normal direction
 *  This is for 2D only
 */
registerMooseObject("farmsApp", ADInitialDamageCycleSim2D);

InputParameters
ADInitialDamageCycleSim2D::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Material used in defining initial damage profile");
  params.addRequiredParam<Real>("len_of_fault", "length of the fault");
  params.addRequiredParam<Real>("sigma", "decay rate");
  params.addRequiredParam<Real>("peak_val", "peak value of the initial damage");
  params.addParam<bool>("use_damage_perturb", false, "Whether to use damage perturbation");
  params.addParam<bool>("use_background_randalpha", false, "Whether to use random alpha background");
  params.addCoupledVar("randalpha", "Random alpha auxiliary variable");
  params.addParam<MaterialPropertyName>("damage_perturb", "",
                                        "Coupled material property for damage perturbation. "
                                        "Must be specified if use_damage_perturb is true.");
  return params;
}

ADInitialDamageCycleSim2D::ADInitialDamageCycleSim2D(const InputParameters & parameters)
  : ADMaterial(parameters),
  _initial_damage(declareADProperty<Real>("initial_damage")),
  _len_of_fault(getParam<Real>("len_of_fault")),
  _sigma(getParam<Real>("sigma")),
  _peak_val(getParam<Real>("peak_val")),
  _use_background_randalpha(getParam<bool>("use_background_randalpha")),
  _randalpha(_use_background_randalpha ? &adCoupledValue("randalpha") : nullptr),
  _use_damage_perturb(getParam<bool>("use_damage_perturb")),
  _damage_perturbation(_use_damage_perturb ? &getADMaterialProperty<Real>("damage_perturb") : nullptr)
{
}

void
ADInitialDamageCycleSim2D::initQpStatefulProperties()
{
}

void
ADInitialDamageCycleSim2D::computeQpProperties()
{

  ADReal x_coord = _q_point[_qp](0); //along the strike direction
  ADReal y_coord = _q_point[_qp](1); //along the normal direction

  ADReal alpha_o = 0;

  ADReal r = 0.0;
  ADReal sigma = _sigma;
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

  if (_use_damage_perturb)
  {
    alpha_o = alpha_o + (*_damage_perturbation)[_qp];
  }

  _initial_damage[_qp] = alpha_o; 

}