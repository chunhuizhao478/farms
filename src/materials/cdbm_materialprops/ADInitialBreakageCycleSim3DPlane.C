//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADInitialBreakageCycleSim3DPlane.h"

/**
 *  Created by Chunhui Zhao, Feb 16th, 2025
 *  A 3D generalization of the initial-damage profile material for a rectangular fault plane.
 *  The plane is defined at y=0 with width `len_of_fault` in the x-direction
 *  and `len_of_fault_dip` in the z-direction.
 *
 *  Exponential decay of damage is computed as a function of distance from
 *  the nearest point on that plane.
 */
registerMooseObject("farmsApp", ADInitialBreakageCycleSim3DPlane);

InputParameters
ADInitialBreakageCycleSim3DPlane::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Material for an initial damage profile in 3D with an "
                             "exponential decay away from a rectangular fault plane in x-z.");
  params.addRequiredParam<Real>("len_of_fault_strike",
                                "Length of the fault in the x-direction");
  params.addRequiredParam<Real>("len_of_fault_dip",
                                "Length of the fault in the z-direction (the 'dip' extent)");
  params.addRequiredParam<Real>("sigma",
                                "Decay rate (controls radial spread of damage)");
  params.addRequiredParam<Real>("peak_val",
                                "Peak value of the initial damage at the plane");
  params.addParam<bool>("use_damage_perturb", false,
                        "Whether to use additional damage perturbation");
  params.addParam<bool>("use_background_randalpha", false,
                        "Whether to use a random alpha background");
  params.addCoupledVar("randalpha", "Random alpha auxiliary variable");
  params.addParam<MaterialPropertyName>(
      "damage_perturb",
      "",
      "Coupled material property for damage perturbation. Must be specified if use_damage_perturb = true.");
  params.addRequiredParam<std::vector<Real>>("nucl_center", "nucleation center (x,y,z)");
  return params;
}

ADInitialBreakageCycleSim3DPlane::ADInitialBreakageCycleSim3DPlane(const InputParameters & parameters)
  : ADMaterial(parameters),
    _initial_breakage(declareADProperty<Real>("initial_breakage")),
    _len_of_fault_strike(getParam<Real>("len_of_fault_strike")),
    _len_of_fault_dip(getParam<Real>("len_of_fault_dip")),
    _sigma(getParam<Real>("sigma")),
    _peak_val(getParam<Real>("peak_val")),
    _use_background_randalpha(getParam<bool>("use_background_randalpha")),
    _randalpha(_use_background_randalpha ? &adCoupledValue("randalpha") : nullptr),
    _use_damage_perturb(getParam<bool>("use_damage_perturb")),
    _damage_perturbation(_use_damage_perturb ? &getADMaterialProperty<Real>("damage_perturb") : nullptr),
    _nucl_center(getParam<std::vector<Real>>("nucl_center"))
{
}

void
ADInitialBreakageCycleSim3DPlane::initQpStatefulProperties()
{
}

void
ADInitialBreakageCycleSim3DPlane::computeQpProperties()
{
  // Coordinates of the current quadrature point
  const ADReal x_coord = _q_point[_qp](0);
  const ADReal y_coord = _q_point[_qp](1);
  const ADReal z_coord = _q_point[_qp](2);

  // Define our rectangular fault plane at y=0.
  // The rectangle is: x in [-0.5*len_of_fault, +0.5*len_of_fault],
  //                  z in [-0.5*len_of_fault_dip+center_point, +0.5*len_of_fault_dip+center_point].

  // 1) Clamp x to the fault rectangle in the x-direction
  const Real x_min = _nucl_center[0] - 0.5 * _len_of_fault_strike;
  const Real x_max = _nucl_center[0] + 0.5 * _len_of_fault_strike;
  const ADReal x_clamped = std::max(x_min, std::min(x_coord, x_max));;

  // 2) Clamp z to the fault rectangle in the z-direction
  const Real z_min = _nucl_center[2] - 0.5 * _len_of_fault_dip;
  const Real z_max = _nucl_center[2] + 0.5 * _len_of_fault_dip;
  const ADReal z_clamped = std::max(z_min, std::min(z_coord, z_max));

  // The plane is at y = 0, so the closest point on the plane to (x, y, z)
  // is (x_clamped, 0, z_clamped). Compute the distance from that point.
  const ADReal dx = x_coord - x_clamped;
  const ADReal dy = y_coord - _nucl_center[1];           // because plane is at y=0
  const ADReal dz = z_coord - z_clamped;
  const ADReal r = std::sqrt(dx * dx + dy * dy + dz * dz);

  // Exponential decay
  ADReal alpha_o = _peak_val * std::exp(-1.0 * (r * r) / (_sigma * _sigma));

  // Optionally include background random alpha
  if (_use_background_randalpha)
    alpha_o = std::max(alpha_o, (*_randalpha)[_qp]);

  // Optionally include additional damage perturbation
  if (_use_damage_perturb)
    alpha_o += (*_damage_perturbation)[_qp];

  // Final initial-breakage value
  _initial_breakage[_qp] = alpha_o;
}