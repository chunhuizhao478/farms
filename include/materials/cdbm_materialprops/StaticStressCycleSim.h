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
 *  Created by Chunhui Zhao, Nov 5th, 2024
 *  Material used in extracting the static stress after the first solve
 */
class StaticStressCycleSim : public Material
{
public:
  static InputParameters validParams();

  StaticStressCycleSim(const InputParameters & parameters);

  virtual void computeQpProperties() override;

  virtual void initQpStatefulProperties() override;

protected:

  /// Material property defining the static stress
  MaterialProperty<RankTwoTensor> & _static_stress;

  /// Material property defining the static stress from the previous time step
  const MaterialProperty<RankTwoTensor> & _static_stress_old;

  /// Material property defining the PK2 stress
  const MaterialProperty<RankTwoTensor> & _pk2_stress;

};