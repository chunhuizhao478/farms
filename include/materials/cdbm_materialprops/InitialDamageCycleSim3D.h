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
 *  Created by Chunhui Zhao, Nov 28th, 2024
 *  Material used in defining initial damage profile with exponential decay along the normal direction
 *  This is for 3D only
 */
class InitialDamageCycleSim3D : public Material
{
public:
  static InputParameters validParams();

  InitialDamageCycleSim3D(const InputParameters & parameters);

  virtual void computeQpProperties() override;

  virtual void initQpStatefulProperties() override;

protected:

  /// Material property initial damage profile
  MaterialProperty<Real> & _initial_damage;

  Real _len_of_fault; // length of the fault
  Real _len_along_dip; // length along fault dip direction
  Real _sigma; // decay rate
  Real _peak_val; // peak value of the initial damage
  /// nucleation center (x,y,z)
  std::vector<Real> _nucl_center;

};