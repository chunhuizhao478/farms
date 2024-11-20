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
#include "FEProblemBase.h"

/**
 *  Created by Chunhui Zhao, Nov 19th, 2024
 *  Material used in get the flag for force damping based on the velocity
 */
class ForceDampingFlag : public Material
{
public:
  static InputParameters validParams();

  ForceDampingFlag(const InputParameters & parameters);

  virtual void computeQpProperties() override;

  virtual void initQpStatefulProperties() override;

protected:

  /// Material property defining the flag for force damping
  MaterialProperty<Real> & _flag;
  const MaterialProperty<Real> & _flag_old;

  /// velocity threshold for force damping
  const Real _vel_maximum_threshold;
  const Real _vel_minimum_threshold;
  const PostprocessorValue & _max_vel_x;  // Add postprocessor value
  const PostprocessorValue & _max_vel_y; 
};