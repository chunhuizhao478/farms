//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ComputeLagrangianStressPK2.h"
#include "GuaranteeConsumer.h"

/// Poroelastic St. Venant-Kirchhoff hyperelasticity
///
/// St. Venant-Kirchhoff hyperelasticity with pore pressure coupling
/// derived from the strain energy function W = lambda / 2 tr(E)^2 + mu tr(E^2)
/// with effective stress approach using Biot's coefficient
///
/// This extends the basic hyperelastic model to include poroelastic effects
/// where the effective stress is computed as:
/// sigma_eff = sigma_total - alpha * p * I
/// where alpha is the Biot coefficient and p is the pore pressure
///
class ComputePoroStVenantKirchhoffStress : public ComputeLagrangianStressPK1, public GuaranteeConsumer
{
public:
  static InputParameters validParams();
  ComputePoroStVenantKirchhoffStress(const InputParameters & parameters);

protected:
  /// Setup function, used to check on isotropy
  virtual void initialSetup() override;

  /// Initial properties
  virtual void initQpStatefulProperties() override;
  
  /// Actual stress/Jacobian update
  virtual void computeQpPK2Stress();

  /// Wrap PK2 -> PK1
  virtual void computeQpPK1Stress() override;

protected:

  /// Elasticity tensor name and property
  const MaterialPropertyName _elasticity_tensor_name;
  const MaterialProperty<RankFourTensor> & _elasticity_tensor;

  /// Pore pressure value
  const VariableValue & _pore_pressure;

  /// Total Lagrange Strain Tensor
  MaterialProperty<RankTwoTensor> & _E;

  /// 2nd PK Stress
  MaterialProperty<RankTwoTensor> & _S;

  /// 2nd PK Tangent (dS/dF)
  MaterialProperty<RankFourTensor> & _C;

  /// off diagnoal jacobian of PK1 with respect to pore pressure
  MaterialProperty<RankTwoTensor> & _pk1_off_diag_jacobian;

  /// First Elastic Strain Invariant
  MaterialProperty<Real> & _I1;

  /// First Elastic Strain Invariant for small strain case
  MaterialProperty<Real> & _I1_small;

  /// derivative of the lagrarangian strain with respect to deformation gradient
  MaterialProperty<RankTwoTensor> & _dI1dF;
  
  /// Material parameters for poroelastic coupling
  const MaterialProperty<Real> & _biot_coeff_eff;
};
