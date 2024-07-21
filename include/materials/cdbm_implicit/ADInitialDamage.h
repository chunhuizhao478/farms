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
 *  Created by Chunhui Zhao, Jul 18th, 2024
 *  ADMaterial used in defining initial damage profile
 */
class ADInitialDamage : public ADMaterial
{
public:
  static InputParameters validParams();

  ADInitialDamage(const InputParameters & parameters);

  virtual void computeQpProperties() override;

  virtual void initQpStatefulProperties() override;

protected:

  /// Material property initial damage profile
  ADMaterialProperty<Real> & _initial_damage;

  /// Material property initial damage profile old
  const MaterialProperty<Real> & _initial_damage_old;

};