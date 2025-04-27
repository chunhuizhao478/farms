//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Material.h"
#include "MooseMesh.h"
#include "MooseTypes.h"

/**
 *  Material used in damage-breakage large deformation formulation, consider full damage evolution equation with diffusion
 *  Created by Chunhui Zhao, Dec 24th, 2024
 */
class DiffusedDamageBreakageMaterialSubApp : public Material
{
public:
  static InputParameters validParams();

  DiffusedDamageBreakageMaterialSubApp(const InputParameters & parameters);

  virtual void initQpStatefulProperties() override;

  virtual void computeQpProperties() override;

  virtual void computegammar(); //compute gamma_r

  virtual void computeStrainRateCd(); //compute strain rate Cd

protected:
  
  //declare material properties:
  MaterialProperty<Real> & _lambda_o_mat; //initial lambda constant value
  MaterialProperty<Real> & _shear_modulus_o_mat; //initial shear modulus value
  MaterialProperty<Real> & _gamma_damaged_r_mat; //material parameter: compliance or fluidity of the fine grain granular material
  MaterialProperty<Real> & _xi_1_mat; //strain invariants ratio: onset of damage evolution
  MaterialProperty<Real> & _xi_0_mat; //strain invariants ratio: onset of damage evolution
  MaterialProperty<Real> & _xi_d_mat; //strain invariants ratio: onset of breakage healing
  MaterialProperty<Real> & _xi_min_mat; //strain invariants ratio: minimum allowable value
  MaterialProperty<Real> & _xi_max_mat; //strain invariants ratio: maximum allowable value  
  MaterialProperty<Real> & _Cd_mat; //coefficient gives positive damage evolution
  MaterialProperty<Real> & _CdCb_multiplier_mat; //multiplier between Cd and Cb
  MaterialProperty<Real> & _CBH_constant_mat; //constant CBH value
  MaterialProperty<Real> & _beta_width_mat; //coefficient gives width of transitional region
  MaterialProperty<Real> & _C1_mat; //coefficient of healing for damage evolution
  MaterialProperty<Real> & _C2_mat; //coefficient of healing for damage evolution
  MaterialProperty<Real> & _D_diffusion_mat; //material parameter: compliance or fluidity of the fine grain granular material
  MaterialProperty<Real> & _initial_damage_mat; //initial damage value
  MaterialProperty<Real> & _I2_mat; //second elastic strain invariant
  MaterialProperty<Real> & _xi_mat; //strain invariants ratio
  MaterialProperty<Real> & _structural_stress_coefficient_mat; //structral_stress_coefficient

  //get const values
  Real _lambda_o_value;
  Real _shear_modulus_o_value;
  Real _xi_0_value;
  Real _xi_d_value;
  Real _xi_min_value;
  Real _xi_max_value;
  Real _Cd_constant_value;
  Real _CdCb_multiplier_value;
  Real _CBH_constant_value;
  Real _beta_width_value;
  Real _C1_value;
  Real _C2_value;
  Real _D_diffusion_value;
  
  //get coupled variables
  const VariableValue & _I2_aux; //second elastic strain invariant
  const VariableValue & _xi_aux; //strain invariants ratio
  const VariableValue & _initial_damage_aux; //initial damage value

  //strain rate dependent Cd options
  bool _use_cd_strain_dependent; //option to use strain rate dependent
  Real _strain_rate_hat; //strain rate for strain-dependent Cd
  Real _cd_hat; //Cd value for strain-dependent Cd
  Real _m_exponent; //exponent for strain-dependent Cd
  const VariableValue & _strain_rate; //strain rate

};