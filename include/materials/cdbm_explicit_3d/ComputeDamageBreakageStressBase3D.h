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
#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include "RotationTensor.h"
#include "ComputeGeneralDamageBreakageStressBase3D.h"

/**
 * ComputeDamageBreakageStressBase3D is the base class for stress tensors
 * computed from MOOSE's strain calculators.
 */
class ComputeDamageBreakageStressBase3D : public ComputeGeneralDamageBreakageStressBase3D
{
public:
  static InputParameters validParams();

  ComputeDamageBreakageStressBase3D(const InputParameters & parameters);
};
