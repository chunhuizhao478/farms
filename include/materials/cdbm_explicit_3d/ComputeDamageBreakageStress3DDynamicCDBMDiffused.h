//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ComputeStressBase.h"

/**
 * ComputeDamageBreakageStress3DDynamicCDBMDiffused put everything inside the computeQpstress without defining
 * additional functions
 
 */
class ComputeDamageBreakageStress3DDynamicCDBMDiffused : public ComputeStressBase
{
public:
  static InputParameters validParams();

  ComputeDamageBreakageStress3DDynamicCDBMDiffused(const InputParameters & parameters);

  virtual void initialSetup() override;

protected:
  virtual void initQpStatefulProperties() override;
  virtual void computeQpStress() override;

  /// @brief Compute gamma_r
  /// @return gamma_r
  Real computegammar();

  /// @brief Compute breakage coefficients
  /// @param gamma_damaged_r
  /// @return a0 a1 a2 a3 in a vector
  std::vector<Real> computecoefficients(Real gamma_damaged_r);

  /// @brief Compute first root of hessian matrix
  /// @param xi 
  /// @return the first root of critical alpha_cr
  Real alphacr_root1(Real xi, Real gamma_damaged_r);

  /// @brief Compute second root of hessian matrix
  /// @param xi 
  /// @return the second root of critical alpha_cr
  Real alphacr_root2(Real xi, Real gamma_damaged_r);

  /// @brief Compute elasticity tensor for small strain
  virtual void computeQpTangentModulus(RankFourTensor & tangent, Real I1, Real I2, Real xi, RankTwoTensor Ee);

  /// @brief Compute the deviatoric strain rate tensor
  virtual void computeDeviatroicStrainRateTensor();

  /// additional variables
  /// initial lambda constant value
  Real _lambda_o;

  /// initial shear modulus value
  Real _shear_modulus_o;

  /// strain invariants ratio: onset of damage evolution
  Real _xi_0;

  /// strain invariants ratio: onset of breakage healing
  Real _xi_d;

  /// energy ratio
  Real _chi;

  /// material parameter: compliance or fluidity of the fine grain granular material
  Real _C_g;

  /// coefficient of power law indexes
  Real _m1;

  /// coefficient of power law indexes
  Real _m2;

  /// alpha_damagedvar_aux
  const VariableValue & _alpha_damagedvar_aux;

  /// B_damagedvar_aux
  const VariableValue & _B_damagedvar_aux;

  /// deviatoric strain rate
  MaterialProperty<Real> & _deviatroic_strain_rate;

  /// deviatroic stress
  MaterialProperty<RankTwoTensor> & _sigma_d;
  const MaterialProperty<RankTwoTensor> & _sigma_d_old;

  /// plastic strain
  MaterialProperty<RankTwoTensor> & _eps_p;
  const MaterialProperty<RankTwoTensor> & _eps_p_old;

  /// total strain
  MaterialProperty<RankTwoTensor> & _eps_total;
  const MaterialProperty<RankTwoTensor> & _eps_total_old;
  
  /// elastic strain
  MaterialProperty<RankTwoTensor> & _eps_e;

  /// I1
  MaterialProperty<Real> & _I1;

  /// I2
  MaterialProperty<Real> & _I2;

  /// xi
  MaterialProperty<Real> & _xi;

  /// lambda
  MaterialProperty<Real> & _lambda;

  /// shear_modulus
  MaterialProperty<Real> & _shear_modulus;

  /// damaged_modulus
  MaterialProperty<Real> & _gamma_damaged;

  /// shear stress perturbation
  const MaterialProperty<Real> & _stress_perturbation;
};