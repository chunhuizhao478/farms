//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ComputeDamageBreakageStressBase3D.h"

/**
 * ComputeDamageBreakageStress3Ddebug put everything inside the computeQpstress without defining
 * additional functions
 
 */
class ComputeDamageBreakageStress3Ddebug : public ComputeDamageBreakageStressBase3D
{
public:
  static InputParameters validParams();

  ComputeDamageBreakageStress3Ddebug(const InputParameters & parameters);

  virtual void initialSetup() override;

protected:
  virtual void initQpStatefulProperties() override;
  virtual void computeQpStress() override;

  /// Function: Compute initial strain based on initial stress
  void setupInitial();

  /// additional variables
  /// strain invariants ratio: onset of damage evolution
  Real _xi_0;

  /// strain invariants ratio: onset of breakage healing
  Real _xi_d;

  /// critical point of three phases
  Real _xi_1;

  /// strain invariants ratio: minimum allowable value
  Real _xi_min;

  /// strain invariants ratio: maximum allowable value
  Real _xi_max;

  /// parameters in granular states
  Real _a0;
  Real _a1;
  Real _a2;
  Real _a3;

  /// coefficient of damage solid modulus
  Real _gamma_damaged_r;

  /// material parameter: compliance or fluidity of the fine grain granular material
  Real _C_g;

  /// coefficient of power law indexes
  Real _m1;

  /// coefficient of power law indexes
  Real _m2;

  /// get old parameters : see definitions in "ComputeGeneralDamageBreakageStressBase"
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
  const MaterialProperty<RankTwoTensor> & _sigma_d_old;

  //add grad term
  const VariableValue & _alpha_grad_x;
  const VariableValue & _alpha_grad_y;
  const VariableValue & _alpha_grad_z;

  /// density
  const MaterialProperty<Real> & _density_old;

  /// diffusion coefficient
  Real _D;

  /// Get initial values
  // const MaterialProperty<RankTwoTensor> & _static_initial_stress_tensor;
  // const MaterialProperty<RankTwoTensor> & _static_initial_strain_tensor;
  // const MaterialProperty<Real> & _I1_initial;
  // const MaterialProperty<Real> & _I2_initial;
  // const MaterialProperty<Real> & _xi_initial;
  const MaterialProperty<Real> & _initial_damage;

  /// damage perturbation
  const MaterialProperty<Real> & _damage_perturbation;

  /// coefficient of positive damage evolution
  Real _Cd_constant;

  /// coefficient of healing of damage evolution
  Real _C1;

  /// coefficient of healing of damage evolution
  Real _C2;

  /// coefficient of width of transitional region
  Real _beta_width;

  /// coefficient of multiplier between Cd and Cb
  Real _CdCb_multiplier;

  /// coefficient of CBH constant
  Real _CBH_constant;

  //initial stress tensor //added
  const MaterialProperty<RankTwoTensor> & _static_initial_stress_tensor;
};