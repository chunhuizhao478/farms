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

/**
 *  Material used in damage-breakage large deformation formulation
 *  Created by Chunhui Zhao, Aug 5th, 2024
 */
class DamageBreakageMaterial : public Material
{
public:
  static InputParameters validParams();

  DamageBreakageMaterial(const InputParameters & parameters);

  virtual void initQpStatefulProperties() override;

  virtual void computeQpProperties() override;

  virtual Real updatedamage(); //update damage variable and return 

  virtual void updatebreakage(); //update breakage variable and return

  virtual void updatemodulus(Real alpha_updated); //update modulus

  virtual void computegammar(); //compute gamma_r

  virtual void computecoefficients(); //compute coefficients a_0 a_1 a_2 a_3 of granular phase

  virtual Real alphacr_root1(Real xi); //alpha cr root 1

  virtual Real alphacr_root2(Real xi); //alpha cr root 2

protected:
  
  //declare material properties:
  MaterialProperty<Real> & _alpha_damagedvar; //damage variable
  MaterialProperty<Real> & _B_breakagevar;    //breakage variable

  MaterialProperty<Real> & _lambda;           //first lamé constant
  MaterialProperty<Real> & _shear_modulus;    //shear modulus
  MaterialProperty<Real> & _damaged_modulus;  //damaged modulus

  MaterialProperty<Real> & _gamma_damaged_r;  //maximum damage modulus

  MaterialProperty<Real> & _a0;               //a0
  MaterialProperty<Real> & _a1;               //a1
  MaterialProperty<Real> & _a2;               //a2
  MaterialProperty<Real> & _a3;               //a3

  MaterialProperty<Real> & _C_g;              //C_g
  MaterialProperty<Real> & _m1;               //m1
  MaterialProperty<Real> & _m2;               //m2
  
  //get old material properties:
  const MaterialProperty<Real> & _alpha_damagedvar_old; //old damage variable
  const MaterialProperty<Real> & _B_breakagevar_old;     //old breakage variable
  const MaterialProperty<Real> & _I2_old;     //old second strain invariants
  const MaterialProperty<Real> & _xi_old;     //old strain invariant ratio
  const MaterialProperty<Real> & _initial_damage; //initial damage 
  const MaterialProperty<Real> & _gamma_damaged_r_old; //old maximum damage modulus
  const MaterialProperty<Real> & _a0_old;     //old a0
  const MaterialProperty<Real> & _a1_old;     //old a1
  const MaterialProperty<Real> & _a2_old;     //old a2
  const MaterialProperty<Real> & _a3_old;     //old a3


  //get const values
  Real _lambda_o;
  Real _shear_modulus_o;
  /// strain invariants ratio: onset of damage evolution
  Real _xi_0;
  /// strain invariants ratio: onset of breakage healing
  Real _xi_d;
  /// strain invariants ratio: minimum allowable value
  Real _xi_min;
  /// strain invariants ratio: maximum allowable value
  Real _xi_max;
  /// coefficient of energy ratio Fb/Fs = chi < 1
  Real _chi;
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
  /// compliance or fluidity of the fine grain granular material
  Real _C_g_value;
  /// coefficient of power law indexes
  Real _m1_value;
  /// coefficient of power law indexes
  Real _m2_value;

};