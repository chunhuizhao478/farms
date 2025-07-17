//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Kernel.h"

class FluidDiffusionGranular : public Kernel
{
public:
  static InputParameters validParams();

  FluidDiffusionGranular(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  /// Helper method to compute gradient of Jp
  RealGradient computeGradJp();

  /// Permeabilities
  const MaterialProperty<Real> & _perm_g;

  /// Fluid viscosity
  const MaterialProperty<Real> & _viscosity;

  /// Deformation gradient
  const MaterialProperty<RankTwoTensor> & _F;

  /// Derivative of Jp with respect to F
  const MaterialProperty<RankTwoTensor> & _dJp_dF;

  /// Plastic Jacobian
  const MaterialProperty<Real> & _Jp;

    /// Derivative of Dp with respect to p
  const MaterialProperty<Real> & _dJp_dp;

  /// Number of nodes in current element
  unsigned int _n_nodes;
};