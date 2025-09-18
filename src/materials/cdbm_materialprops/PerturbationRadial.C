//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PerturbationRadial.h"

/**
 *  Created by Chunhui Zhao, Nov 26th, 2024
 *  Material used in Create Time Dependent Damage/Shear Stress Perturbation in the Dynamic Solve
 */
registerMooseObject("farmsApp", PerturbationRadial);

InputParameters
PerturbationRadial::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Material used in Create Time Dependent Damage Perturbation in the Dynamic Solve");
  params.addRequiredParam<std::vector<Real>>("nucl_center", "nucleation center (x,y,z)");
  params.addRequiredParam<Real>("peak_value","the peak damage value apply region normal to the fault (exponential decay)");
  params.addRequiredParam<Real>("thickness","the thickness used in apply damage value normal to the fault (exponential decay)");
  params.addRequiredParam<Real>("length","the standard deviation used in apply damage value normal to the fault (exponential decay)");
  params.addRequiredParam<Real>("duration","duration to reach peak damage");
  params.addRequiredParam<std::string>("perturbation_type", "Type of perturbation: 'damage' or 'shear_stress'"); // New parameter
  params.addParam<Real>("sigma_divisor", 2.0, "sigma value = (length / sigma_divisor)");
  return params;
}

PerturbationRadial::PerturbationRadial(const InputParameters & parameters)
  : Material(parameters),
  _damage_perturbation(declareProperty<Real>("damage_perturbation")),
  _damage_perturbation_old(getMaterialPropertyOldByName<Real>("damage_perturbation")),
  _shear_stress_perturbation(declareProperty<Real>("shear_stress_perturbation")),
  _shear_stress_perturbation_old(getMaterialPropertyOldByName<Real>("shear_stress_perturbation")),
  _nucl_center_mat(declareProperty<std::vector<Real>>("nucl_center_mat")),
  _thickness_mat(declareProperty<Real>("thickness_mat")),
  _length_mat(declareProperty<Real>("length_mat")),
  _nucl_center(getParam<std::vector<Real>>("nucl_center")),
  _peak_value(getParam<Real>("peak_value")),
  _thickness(getParam<Real>("thickness")),
  _length(getParam<Real>("length")),
  _duration(getParam<Real>("duration")),
  _perturbation_type(getParam<std::string>("perturbation_type")), // Initialize new parameter
  _sigma_divisor(getParam<Real>("sigma_divisor"))
{
  //in case I'm stupid
  if (_nucl_center.size() != 3){
    mooseError("fault plane parameter must accepts 3 numbers!");
  }
}

void
PerturbationRadial::initQpStatefulProperties()
{
  _damage_perturbation[_qp] = 0.0;
  _shear_stress_perturbation[_qp] = 0.0;
  _nucl_center_mat[_qp] = _nucl_center;
  _thickness_mat[_qp] = _thickness;
  _length_mat[_qp] = _length;
}

void
PerturbationRadial::computeQpProperties()
{

  // Get the current point coordinates in the mesh
  const Real xcoord = _q_point[_qp](0); // strike direction
  const Real ycoord = _q_point[_qp](1); // normal direction
  const Real zcoord = _q_point[_qp](2); // dip direction

  // We define a 2D Gaussian in the XZ plane, ignoring y in the exponent
  // The characteristic "sigma" is length / sigma_divisor
  const Real sigma_x = _length / _sigma_divisor;
  const Real sigma_z = _length / _sigma_divisor;
  const Real gaussian_factor = _peak_value; // maximum amplitude of the Gaussian

  // Distance in XZ from the nucleation center
  // If your center is (0, 0, -7500), then _nucl_center might be [0, 0, -7500].
  const Real dx = xcoord - _nucl_center[0];
  const Real dz = zcoord - _nucl_center[2];

  // 2D Gaussian distribution in the XZ plane:
  //    G(x,z) = peak_value * exp( - [dx^2 / (2*sigma_x^2) + dz^2 / (2*sigma_z^2)] )
  const Real gaussian_value = gaussian_factor *
                             std::exp(-((dx * dx) / (2.0 * sigma_x * sigma_x) +
                                        (dz * dz) / (2.0 * sigma_z * sigma_z)));

  // Scale Gaussian value over time
  Real scaled_gaussian_value = 0.0;
  if (_t <= _duration)
  {
    scaled_gaussian_value = gaussian_value * (_t / _duration);
  }
  else
  {
    scaled_gaussian_value = gaussian_value;
  }

  // Check Z direction constraint and apply perturbation
  Real dalpha_damage = 0.0;
  Real dalpha_stress = 0.0;
  if (ycoord >= _nucl_center[1] - _thickness / 2.0 && ycoord <= _nucl_center[1] + _thickness / 2.0)
  {
    dalpha_damage = scaled_gaussian_value;
    dalpha_stress = scaled_gaussian_value;
  }
  else
  {
    dalpha_damage = _damage_perturbation_old[_qp];
    dalpha_stress = _shear_stress_perturbation_old[_qp];
  }

  if (_perturbation_type == "damage")
  {
    _damage_perturbation[_qp] = dalpha_damage;
    _shear_stress_perturbation[_qp] = 0.0;
  }
  else if (_perturbation_type == "shear_stress")
  {
    _damage_perturbation[_qp] = 0.0;
    _shear_stress_perturbation[_qp] = dalpha_stress;
  }
  else
  {
    mooseError("Invalid perturbation type: " + _perturbation_type);
  }

  // Update nucleation center, thickness, and length
  _nucl_center_mat[_qp] = _nucl_center;
  _thickness_mat[_qp] = _thickness;
  _length_mat[_qp] = _length;

}