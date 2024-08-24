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
 *  Created by Chunhui Zhao, Aug 14th, 2024
 *  Material used in defining initial damage profile with exponential decay along the normal direction
 */
class InitialDamageBenchmarkSperical : public Material
{
public:
  static InputParameters validParams();

  InitialDamageBenchmarkSperical(const InputParameters & parameters);

  virtual void computeQpProperties() override;

  virtual void initQpStatefulProperties() override;

protected:

  /// Material property initial damage profile
  MaterialProperty<Real> & _initial_damage;

  /// nucleation center (x,y,z)
  std::vector<Real> _nucl_center;

  /// peak damage (exponential decay)
  Real _peak_damage;

  /// sigma (exponential decay)
  Real _sigma;

};