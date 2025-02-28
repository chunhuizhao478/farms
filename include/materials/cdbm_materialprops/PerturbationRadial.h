//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Material.h"

/**
 *  Created by Chunhui Zhao, Nov 26th, 2024
 *  Material used in Create Time Dependent Damage/Shear Stress Perturbation in the Dynamic Solve
 */
class PerturbationRadial : public Material
{
public:
  static InputParameters validParams();

  PerturbationRadial(const InputParameters & parameters);

  virtual void computeQpProperties() override;

  virtual void initQpStatefulProperties() override;

protected:

  /// Material property initial damage profile
  MaterialProperty<Real> & _damage_perturbation;

  /// Material property old initial damage profile
  const MaterialProperty<Real> & _damage_perturbation_old;

  /// Material property breakage perturbation
  MaterialProperty<Real> & _breakage_perturbation;

  /// Material property old breakage perturbation
  const MaterialProperty<Real> & _breakage_perturbation_old;

  /// Material property shear stress perturbation
  MaterialProperty<Real> & _shear_stress_perturbation;

  /// Material property old shear stress perturbation
  const MaterialProperty<Real> & _shear_stress_perturbation_old;

  /// Material property nucleation center
  MaterialProperty<std::vector<Real>> & _nucl_center_mat;

  /// Material property thickness
  MaterialProperty<Real> & _thickness_mat;

  /// Material property length
  MaterialProperty<Real> & _length_mat;

  /// nucleation center (x,y,z)
  std::vector<Real> _nucl_center;

  /// peak damage (exponential decay)
  Real _peak_value;

  /// thickness
  Real _thickness;

  /// length
  Real _length;

  /// duration to reach peak damage
  Real _duration;

  std::string _perturbation_type; // New parameter to choose perturbation type

  /// sigma value (length / coefficient)
  Real _sigma_divisor;

};