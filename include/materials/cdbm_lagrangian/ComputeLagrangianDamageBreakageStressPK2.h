//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ComputeLagrangianStressPK1.h"

/// Native interface for providing the 2nd Piola Kirchhoff stress
///
/// This class *implements* the 2nd PK stress update, providing:
///   1) The 2nd PK stress
///   2) The derivative of the 2nd PK stress wrt the Cauchy-Green strain
///
/// and wraps these to provide:
///   1) The 1st PK stress
///   2) d(PK1)/d(F)
///
class ComputeLagrangianDamageBreakageStressPK2 : public ComputeLagrangianStressPK1
{
public:
  static InputParameters validParams();
  ComputeLagrangianDamageBreakageStressPK2(const InputParameters & parameters);

protected:
  /// Initial properties
  virtual void initQpStatefulProperties() override;
  /// Wrap PK2 -> PK1
  virtual void computeQpPK1Stress() override;
  /// Provide the PK2 stress and dPK2/dC
  virtual void computeQpPK2Stress();
  /// @brief Update plastic deformation gradient
  /// @return plastic deformation gradient
  virtual RankTwoTensor computeQpFp();
  /// @brief Compute tangent modulus components vector
  /// @return tangent modulus components vector, we pass reference and modify in-place
  virtual void computeQpTangentModulus(RankFourTensor & tangent, Real I1, Real I2, Real xi, RankTwoTensor Ee);

protected:
  /* Declare Material Properties */
  /// Plastic Deformation Gradient
  MaterialProperty<RankTwoTensor> & _Fp;
  /// Determinant of Plastic Deformation Gradient
  MaterialProperty<Real> & _Jp;
  /// Elastic Deformation Gradient
  MaterialProperty<RankTwoTensor> & _Fe;
  /// Deviatoric Stress Tensor
  MaterialProperty<RankTwoTensor> & _Tau;
  /// Green-Lagrange Elastic Strain Tensor
  MaterialProperty<RankTwoTensor> & _Ee;
  /// Plastic Strain Tensor
  MaterialProperty<RankTwoTensor> & _Ep;
  /// Total Lagrange Strain Tensor
  MaterialProperty<RankTwoTensor> & _E;
  /// First Elastic Strain Invariant
  MaterialProperty<Real> & _I1;
  /// Second Elastic Strain Invariant
  MaterialProperty<Real> & _I2;
  /// Strain Invariant Ratio
  MaterialProperty<Real> & _xi;
  /// 2nd PK Stress
  MaterialProperty<RankTwoTensor> & _S;
  /// 2nd PK Tangent (dS/dF)
  MaterialProperty<RankFourTensor> & _C;

  /* Get Up-to-date Material Properties*/
  /// Lambda
  const MaterialProperty<Real> & _lambda_const;
  /// Shear Modulus
  const MaterialProperty<Real> & _shear_modulus;
  /// Damage Modulus
  const MaterialProperty<Real> & _damaged_modulus;
  /// Breakage Variable
  const MaterialProperty<Real> & _B_breakagevar;
  
  /* Get Old Material Properties */
  /// Breakage Variable
  const MaterialProperty<Real> & _B_breakagevar_old;
  /// Deviatroic Stress Tensor
  const MaterialProperty<RankTwoTensor> & _Tau_old;
  /// Plastic Deformation Gradient
  const MaterialProperty<RankTwoTensor> & _Fp_old;
  /// Deformation Gradient
  const MaterialProperty<RankTwoTensor> & _F_old;
  /// Plastic Strain
  const MaterialProperty<RankTwoTensor> & _Ep_old;
  
  /* Get Constant Parameters */
  /// material parameter: compliance or fluidity of the fine grain granular material
  const MaterialProperty<Real> & _C_g; 
  /// coefficient of power law indexes
  const MaterialProperty<Real> & _m1;
  /// coefficient of power law indexes
  const MaterialProperty<Real> & _m2;
  /// parameters in granular states
  const MaterialProperty<Real> & _a0;
  const MaterialProperty<Real> & _a1;
  const MaterialProperty<Real> & _a2;
  const MaterialProperty<Real> & _a3;  
  /// dimension of the problem
  const unsigned int _dim;
};