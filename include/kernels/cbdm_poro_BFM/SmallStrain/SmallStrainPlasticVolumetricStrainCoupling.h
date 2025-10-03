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

class SmallStrainPlasticVolumetricStrainCoupling : public Kernel
{
public:
  static InputParameters validParams();

  SmallStrainPlasticVolumetricStrainCoupling(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  /// Number of displacement components
  const unsigned int _ndisp;

  /// Variable numbers for displacement components
  std::vector<unsigned int> _disp_var;

  /// Plastic strain
  const MaterialProperty<RankTwoTensor> & _eps_p;

  /// Plastic strain old
  const MaterialProperty<RankTwoTensor> & _eps_p_old;

  /// Derivatives for off-diagonal Jacobian
  const MaterialProperty<RankTwoTensor> & _deps_p_dp;
  const MaterialProperty<RankFourTensor> & _deps_p_deps;
};

