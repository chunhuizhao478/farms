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
 * ADComputeDamageBreakageStress computes the stress following linear elasticity theory (small
 * strains)
 */
//template <typename RankTwoTensor, typename RankFourTensor, typename Real>
class ADComputeDamageBreakageStress : public ADComputeDamageBreakageStressBase
{
public:
  static InputParameters validParams();

  ADComputeDamageBreakageStress(const InputParameters & parameters);

  virtual void initialSetup() override;

protected:
  virtual void computeQpStress() override;

  //additional parameters

  ADReal _lambda_o;
  ADReal _shear_modulus_o;

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

  /// coefficient of positive damage evolution
  ADReal _Cd_constant;

  /// coefficient of healing of damage evolution
  ADReal _C1;

  /// coefficient of healing of damage evolution
  ADReal _C2;

  /// coefficient of width of transitional region
  ADReal _beta_width;

  /// coefficient of multiplier between Cd and Cb
  ADReal _CdCb_multiplier;

  /// coefficient of CBH constant
  ADReal _CBH_constant;

  /// get material prop
  const MaterialProperty<Real> & _initial_damage;

  /// get old parameters : see definitions in "ADComputeDamageBreakageStressBase"
  const MaterialProperty<Real> & _alpha_damagedvar_old;
  const MaterialProperty<Real> & _B_breakagevar_old;
  const MaterialProperty<Real> & _xi_old;
  const MaterialProperty<Real> & _I2_old;
  const MaterialProperty<Real> & _shear_modulus_old;
  const MaterialProperty<Real> & _gamma_damaged_old;
  const MaterialProperty<RankTwoTensor> & _eps_p_old;
  const MaterialProperty<RankTwoTensor> & _eps_e_old;
  const MaterialProperty<RankTwoTensor> & _sigma_d_old;

};