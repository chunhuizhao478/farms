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
 *  Created by Chunhui Zhao, Jul 17th, 2024
 *  Material used to compute xi of elastic strain
 */
class ComputeXi : public Material
{
public:
  static InputParameters validParams();

  ComputeXi(const InputParameters & parameters);

  virtual void computeQpProperties() override;

protected:

  /// @brief Define the xi as material property
  MaterialProperty<Real> & _xi;

  MaterialProperty<Real> & _I1;

  MaterialProperty<Real> & _I2;

  /// elastic strain
  const MaterialProperty<RankTwoTensor> & _mechanical_strain;

};