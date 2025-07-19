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

  /// Permeabilities
  const MaterialProperty<Real> & _perm_g;

  /// Deformation gradient
  const MaterialProperty<RankTwoTensor> & _F;

  /// Plastic Jacobian
  const MaterialProperty<Real> & _Jp;

    /// Derivative of Dp with respect to p
  const MaterialProperty<Real> & _dJp_dp;
};