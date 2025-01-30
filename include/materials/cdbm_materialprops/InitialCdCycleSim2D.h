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
 *  Created by Chunhui Zhao, Jan 29th, 2025
 *  ADMaterial used in defining initial Cd profile with exponential decay along the normal direction
 *  This is for 2D only
 */
class InitialCdCycleSim2D : public Material
{
public:
  static InputParameters validParams();

  InitialCdCycleSim2D(const InputParameters & parameters);

  virtual void computeQpProperties() override;

  virtual void initQpStatefulProperties() override;

protected:

  /// Material property initial damage profile
  MaterialProperty<Real> & _initial_damage;

  Real _len_of_fault; // length of the fault
  Real _sigma; // decay rate
  Real _peak_val; // peak value of the initial damage

  // Add new members
  const bool _use_background_randalpha;
  const VariableValue * _randalpha;

};