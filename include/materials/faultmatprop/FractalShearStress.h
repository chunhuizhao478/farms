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
 *  Created by Chunhui Zhao, Apr 12th, 2025
 *  Material used in Generate Fractal Shear Stress Distribution
 */
class FractalShearStress : public Material
{
public:
  static InputParameters validParams();

  FractalShearStress(const InputParameters & parameters);

  virtual void computeQpProperties() override;

  virtual void initQpStatefulProperties() override;

protected:

  /// Material property initial damage profile
  MaterialProperty<Real> & _fractal_shear_stress;
  const MaterialProperty<Real> & _fractal_shear_stress_old;
  const std::string _csv_file;
};