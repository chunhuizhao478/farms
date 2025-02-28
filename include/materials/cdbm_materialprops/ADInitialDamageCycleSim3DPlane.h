//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ADMaterial.h"

/**
 *  Created by Chunhui Zhao, Feb 16th, 2025
 *  A 3D generalization of the initial-damage profile material for a rectangular fault plane.
 *  The plane is defined at y=0 with width `len_of_fault` in the x-direction
 *  and `len_of_fault_dip` in the z-direction.
 *
 *  Exponential decay of damage is computed as a function of distance from
 *  the nearest point on that plane.
 */
class ADInitialDamageCycleSim3DPlane : public ADMaterial
{
public:
  static InputParameters validParams();

  ADInitialDamageCycleSim3DPlane(const InputParameters & parameters);

  virtual void computeQpProperties() override;

  virtual void initQpStatefulProperties() override;

protected:

  /// Material property initial damage profile
  ADMaterialProperty<Real> & _initial_damage;

  /// Length of the fault in the x-direction
  Real _len_of_fault_strike;

  /// Length of the fault in the z-direction (the 'dip' extent)
  Real _len_of_fault_dip;

  /// sigma (exponential decay)
  Real _sigma;

  /// Peak value of the initial damage at the plane
  Real _peak_val;

  /// Whether to use a random alpha background
  const bool _use_background_randalpha;

  /// Random alpha auxiliary variable
  const ADVariableValue * _randalpha;

  /// Whether to use additional damage perturbation
  const bool _use_damage_perturb;

  /// Coupled material property for damage perturbation
  const ADMaterialProperty<Real> * _damage_perturbation;

  /// nucleation center (x,y,z)
  std::vector<Real> _nucl_center;

};