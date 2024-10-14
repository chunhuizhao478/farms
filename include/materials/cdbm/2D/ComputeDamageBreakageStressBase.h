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
#include "ComputeGeneralDamageBreakageStressBase.h"

/**
 * ComputeDamageBreakageStressBase is the base class for stress tensors
 * computed from MOOSE's strain calculators.
 */
class ComputeDamageBreakageStressBase : public ComputeGeneralDamageBreakageStressBase
{
public:
  static InputParameters validParams();

  ComputeDamageBreakageStressBase(const InputParameters & parameters);
};
