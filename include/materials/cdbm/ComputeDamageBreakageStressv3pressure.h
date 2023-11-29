//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ComputeDamageBreakageStressBase.h"

/**
 * ComputeDamageBreakageStressv3pressure put everything inside the computeQpstress without defining
 * additional functions
 
 */
class ComputeDamageBreakageStressv3pressure : public ComputeDamageBreakageStressBase
{
public:
  static InputParameters validParams();

  ComputeDamageBreakageStressv3pressure(const InputParameters & parameters);

  virtual void initialSetup() override;

protected:
  virtual void computeQpStress() override;

  /// Function: Compute initial strain based on initial stress
  void setupInitial();

  /// Name of the elasticity tensor material property
  //const std::string _elasticity_tensor_name;
  /// Elasticity tensor material property
  //const MaterialProperty<RankFourTensor> & _elasticity_tensor;

  //initial stress tensor
  const MaterialProperty<RankTwoTensor> & _static_initial_stress_tensor;

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
  const MaterialProperty<Real> & _eqv_plastic_strain_old; //old eqv plastic strain

  /// updated damage and breakage parameters computed from subApp
  //Note: pass reference(&) instead of value, otherwise it may occur segmentation fault 11 error
  const VariableValue & _alpha_in;
  const VariableValue & _B_in;

  //add grad term
  const VariableValue & _alpha_grad_x;
  const VariableValue & _alpha_grad_y;

  // pore pressure
  Real _effec_sts_coeff;
  const VariableValue &  _pressure; //constant pointer to _pressure //the value being pointed at can change but the pointer can't

  /// density
  const MaterialProperty<Real> & _density;

  /// diffusion coefficient
  Real _D;

  /// Function: deltaij
  Real deltaij(int i, int j);

  /// Function: epsilonij - take component of elastic strain
  Real epsilonij(int i, 
                 int j,
                 Real eps11e_in,
                 Real eps22e_in,
                 Real eps12e_in);

  Real grad_alpha(int i, 
                  Real alpha_grad_x,
                  Real alpha_grad_y);

  /// Function: compute stress components
  Real computeStressComps(int i, 
                          int j,
                          Real xi_in,
                          Real I1_in,
                          Real B_in,
                          Real lambda_in,
                          Real gamma_damaged_in,
                          Real shear_modulus_in,
                          Real eps11e_in,
                          Real eps22e_in,
                          Real eps12e_in,
                          Real alpha_grad_x,
                          Real alpha_grad_y,
                          Real D);
};