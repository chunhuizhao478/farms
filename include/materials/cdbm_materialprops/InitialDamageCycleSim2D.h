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
 *  Created by Chunhui Zhao, Oct 20th, 2024
 *  ADMaterial used in defining initial damage profile with exponential decay along the normal direction
 *  This is for 2D only
 */
class InitialDamageCycleSim2D : public Material
{
public:
  static InputParameters validParams();

  InitialDamageCycleSim2D(const InputParameters & parameters);

  virtual void computeQpProperties() override;

  virtual void initQpStatefulProperties() override;

protected:

  /// Material property initial damage profile
  MaterialProperty<Real> & _initial_damage;

};