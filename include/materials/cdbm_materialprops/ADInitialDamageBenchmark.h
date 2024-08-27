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
 *  Created by Chunhui Zhao, Aug 14th, 2024
 *  ADMaterial used in defining initial damage profile with exponential decay along the normal direction
 */
class ADInitialDamageBenchmark : public ADMaterial
{
public:
  static InputParameters validParams();

  ADInitialDamageBenchmark(const InputParameters & parameters);

  virtual void computeQpProperties() override;

  virtual void initQpStatefulProperties() override;

protected:

  /// Material property initial damage profile
  ADMaterialProperty<Real> & _initial_damage;

  /// nucleation center (x,y,z)
  std::vector<Real> _nucl_center;

  /// region of fault plane (xmin,xmax,ymin,ymax)
  std::vector<Real> _fault_plane;

  /// nucleation distance
  Real _nucl_distance;

  /// nucleation thickness
  Real _nucl_thickness;

  /// nucleation damage
  Real _nucl_damage;

  /// peak damage (exponential decay)
  Real _peak_damage;

  /// sigma (exponential decay)
  Real _sigma;

};