//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ADComputeDamageBreakageStressBase.h"

/**
 * ADComputeDamageBreakageStressDev computes the stress following linear elasticity theory (small
 * strains)
 */
//template <typename RankTwoTensor, typename RankFourTensor, typename Real>
class ADComputeDamageBreakageStressDev : public ADComputeDamageBreakageStressBase
{
public:
  static InputParameters validParams();

  ADComputeDamageBreakageStressDev(const InputParameters & parameters);

  virtual void initialSetup() override;

protected:
  virtual void computeQpStress() override;

  /// Function: Compute initial strain based on initial stress
  void setupInitial();

  /// Name of the elasticity tensor material property
  //const std::string _elasticity_tensor_name;

  /// Elasticity tensor material property
  //const ADMaterialProperty<RankFourTensor> & _elasticity_tensor;

  //additional parameters
  //initial stress tensor
  const ADMaterialProperty<RankTwoTensor> & _static_initial_stress_tensor;

  //initial strain tensor
  const ADMaterialProperty<RankTwoTensor> & _static_initial_strain_tensor;  

  /// additional variables
  /// strain invariants ratio: onset of damage evolution
  ADReal _xi_0;

  /// strain invariants ratio: onset of breakage healing
  ADReal _xi_d;

  /// critical point of three phases
  ADReal _xi_1;

  /// strain invariants ratio: minimum allowable value
  ADReal _xi_min;

  /// strain invariants ratio: maximum allowable value
  ADReal _xi_max;

  /// parameters in granular states
  ADReal _a0;
  ADReal _a1;
  ADReal _a2;
  ADReal _a3;

  /// coefficient of damage solid modulus
  ADReal _gamma_damaged_r;

  /// material parameter: compliance or fluidity of the fine grain granular material
  ADReal _C_g;

  /// coefficient of power law indexes
  ADReal _m1;

  /// coefficient of power law indexes
  ADReal _m2;

  /// get old parameters : see definitions in "adComputeDamageBreakageStressBase"
  const MaterialProperty<Real> & _alpha_damagedvar_old;
  const MaterialProperty<Real> & _B_old;
  const MaterialProperty<Real> & _xi_old;
  const MaterialProperty<Real> & _I1_old;
  const MaterialProperty<Real> & _I2_old;
  const MaterialProperty<Real> & _lambda_old;
  const MaterialProperty<Real> & _shear_modulus_old;
  const MaterialProperty<Real> & _gamma_damaged_old;
  const MaterialProperty<RankTwoTensor> & _eps_total_old;
  const MaterialProperty<RankTwoTensor> & _mechanical_strain_old;
  const MaterialProperty<RankTwoTensor> & _eps_p_old;
  const MaterialProperty<RankTwoTensor> & _eps_e_old;

  /// updated damage and breakage parameters computed from subApp
  //Note: pass reference(&) instead of value, otherwise it may occur segmentation fault 11 error
  const ADVariableValue & _alpha_in;
  const ADVariableValue & _B_in;

  //add grad term
  const ADVariableValue & _alpha_grad_x;
  const ADVariableValue & _alpha_grad_y;
  const ADVariableValue & _alpha_grad_z;

  /// density
  const ADMaterialProperty<Real> & _density_old;

  /// diffusion coefficient
  ADReal _D;
  
  /// initial damage
  const ADVariableValue & _initial_alpha; 
  // usingComputeDamageBreakageStressBaseMembers;
};

//typedef ADComputeDamageBreakageStressDev<RankTwoTensor, RankFourTensor, Real>
//    ADComputeDamageBreakageStressDev;
//typedef ADComputeDamageBreakageStressDev<SymmetricRankTwoTensor, SymmetricRankFourTensor>
//    ADSymmetricDamageBreakageStress;