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
 *  Created by Chunhui Zhao, Jul 17th, 2024
 *  ADMaterial used to compute xi of elastic strain
 */
class ADComputeXi : public ADMaterial
{
public:
  static InputParameters validParams();

  ADComputeXi(const InputParameters & parameters);

  virtual void computeQpProperties() override;

protected:

  /// @brief Define the xi as material property
  ADMaterialProperty<Real> & _xi;

  ADMaterialProperty<Real> & _I1;

  ADMaterialProperty<Real> & _I2;

  /// elastic strain
  const ADMaterialProperty<RankTwoTensor> & _mechanical_strain;

};