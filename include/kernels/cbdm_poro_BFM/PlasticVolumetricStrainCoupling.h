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

class PlasticVolumetricStrainCoupling : public Kernel
{
public:
  static InputParameters validParams();

  PlasticVolumetricStrainCoupling(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  /// Number of displacement components
  const unsigned int _ndisp;

  /// Variable numbers for displacement components
  std::vector<unsigned int> _disp_var;

  /// Plastic Jacobian
  const MaterialProperty<Real> & _Jp;

  /// Plastic strain rate tensor
  const MaterialProperty<RankTwoTensor> & _Dp;

  /// Derivatives for off-diagonal Jacobian
  const MaterialProperty<RankTwoTensor> & _dJp_dF;
  const MaterialProperty<Real> & _dJp_dp;
  const MaterialProperty<RankFourTensor> & _dDp_dF;
  const MaterialProperty<RankTwoTensor> & _dDp_dp;
};