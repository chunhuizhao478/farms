//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "FarmsComputeLagrangianStrainBase.h"

template <>
inline InputParameters
FarmsComputeLagrangianStrainBase<GradientOperatorCartesian>::validParams()
{
  InputParameters params = FarmsComputeLagrangianStrainBase::baseParams();
  params.addClassDescription("Compute strain in Cartesian coordinates.");
  return params;
}

template <>
inline void
FarmsComputeLagrangianStrainBase<GradientOperatorCartesian>::initialSetup()
{
  if (getBlockCoordSystem() != Moose::COORD_XYZ)
    mooseError("This kernel should only act in Cartesian coordinates.");
}

typedef FarmsComputeLagrangianStrainBase<GradientOperatorCartesian> FarmsComputeLagrangianStrain;