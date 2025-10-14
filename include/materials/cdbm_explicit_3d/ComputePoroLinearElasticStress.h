//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ComputeStressBase.h"

/**
 * ComputePoroLinearElasticStress computes the stress following linear elasticity theory (small strains)
 */
class ComputePoroLinearElasticStress : public ComputeStressBase
{
public:
  static InputParameters validParams();

  ComputePoroLinearElasticStress(const InputParameters & parameters);

  virtual void initialSetup() override;

protected:
  virtual void computeQpStress() override;

  /// Name of the elasticity tensor material property
  const std::string _elasticity_tensor_name;
  /// Elasticity tensor material property
  const MaterialProperty<RankFourTensor> & _elasticity_tensor;

  /// pore pressure coupled variable
  const VariableValue & _pore_pressure;

  /// First Elastic Strain Invariant
  MaterialProperty<Real> & _I1;

  /// stress off-diagonal Jacobian
  MaterialProperty<RankTwoTensor> & _stress_off_diag_jacobian;

  /// Material parameters for poroelastic coupling
  const MaterialProperty<Real> & _biot_coeff_eff;
};