//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "LagrangianStressDivergenceBase.h"
#include "GradientOperator.h"

/// Enforce equilibrium with a total Lagrangian formulation
///
/// This class enforces equilibrium when used in conjunction with
/// the corresponding strain calculator (CalculateStrainLagrangianKernel)
/// and with either a stress calculator that provides the
/// 1st PK stress ("pk1_stress") and the derivative of the 1st PK stress
/// with respect to the deformation gradient ("pk1_jacobian")
///
/// This kernel should be used with the new "ComputeLagrangianStressBase"
/// stress update system and the "ComputeLagrangianStrain" system for strains.
///
/// Optional pore pressure coupling is available when include_fluid=true,
/// which adds off-diagonal coupling terms for poromechanical problems.
///
template <class G>
class TotalLagrangianTotalStressDivergenceBase : public LagrangianStressDivergenceBase, G
{
public:
  static InputParameters baseParams();
  static InputParameters validParams();
  TotalLagrangianTotalStressDivergenceBase(const InputParameters & parameters);
  virtual void initialSetup() override;

protected:
  virtual RankTwoTensor gradTest(unsigned int component) override;
  virtual RankTwoTensor gradTrial(unsigned int component) override;
  virtual void precalculateJacobianDisplacement(unsigned int component) override;
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobianDisplacement(unsigned int alpha, unsigned int beta) override;
  /// Compute off-diagonal Jacobian with respect to pore pressure
  virtual Real computeQpJacobianPorePressure(unsigned int jvar);
  virtual Real computeQpJacobianTemperature(unsigned int cvar) override;
  virtual Real computeQpJacobianOutOfPlaneStrain() override;
  
  /// The 1st Piola-Kirchhoff stress
  const MaterialProperty<RankTwoTensor> & _pk1;
  
  /// The derivative of the PK1 stress with respect to the deformation gradient
  const MaterialProperty<RankFourTensor> & _dpk1;

  /// Pore pressure variable (following temperature pattern)
  const VariableValue & _pore_pressure;

  /// Pore pressure variable 
  const unsigned int _pore_pressure_var;
  
  /// PK1 off-diagonal Jacobian with respect to pore pressure
  const MaterialProperty<RankTwoTensor> & _pk1_off_diag_jacobian;

private:
  /// The unstabilized trial function gradient
  virtual RankTwoTensor gradTrialUnstabilized(unsigned int component);
  
  /// The stabilized trial function gradient
  virtual RankTwoTensor gradTrialStabilized(unsigned int component);
};