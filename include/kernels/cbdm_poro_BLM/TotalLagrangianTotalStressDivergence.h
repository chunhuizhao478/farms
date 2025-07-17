//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once
#include "TotalLagrangianTotalStressDivergenceBase.h"

// Cartesian coordinates
template <>
inline InputParameters
TotalLagrangianTotalStressDivergenceBase<GradientOperatorCartesian>::validParams()
{
  InputParameters params = TotalLagrangianTotalStressDivergenceBase::baseParams();
  params.addClassDescription(
      "Enforce equilibrium with a total Lagrangian formulation in Cartesian coordinates with optional pore pressure coupling.");
  return params;
}

template <>
inline void
TotalLagrangianTotalStressDivergenceBase<GradientOperatorCartesian>::initialSetup()
{
  if (getBlockCoordSystem() != Moose::COORD_XYZ)
    mooseError("This kernel should only act in Cartesian coordinates.");
}

// Axisymmetric Cylindrical coordinates  
template <>
inline InputParameters
TotalLagrangianTotalStressDivergenceBase<GradientOperatorAxisymmetricCylindrical>::validParams()
{
  InputParameters params = TotalLagrangianTotalStressDivergenceBase::baseParams();
  params.addClassDescription(
      "Enforce equilibrium with a total Lagrangian formulation in axisymmetric cylindrical coordinates with optional pore pressure coupling.");
  return params;
}

template <>
inline void
TotalLagrangianTotalStressDivergenceBase<GradientOperatorAxisymmetricCylindrical>::initialSetup()
{
  if (getBlockCoordSystem() != Moose::COORD_RZ)
    mooseError("This kernel should only act in axisymmetric cylindrical coordinates.");
}

// Centrosymmetric Spherical coordinates
template <>
inline InputParameters
TotalLagrangianTotalStressDivergenceBase<GradientOperatorCentrosymmetricSpherical>::validParams()
{
  InputParameters params = TotalLagrangianTotalStressDivergenceBase::baseParams();
  params.addClassDescription(
      "Enforce equilibrium with a total Lagrangian formulation in centrosymmetric spherical coordinates with optional pore pressure coupling.");
  return params;
}

template <>
inline void
TotalLagrangianTotalStressDivergenceBase<GradientOperatorCentrosymmetricSpherical>::initialSetup()
{
  if (getBlockCoordSystem() != Moose::COORD_RSPHERICAL)
    mooseError("This kernel should only act in centrosymmetric spherical coordinates.");
}

// Create concrete class name using typedef
typedef TotalLagrangianTotalStressDivergenceBase<GradientOperatorCartesian>
    TotalLagrangianTotalStressDivergence;