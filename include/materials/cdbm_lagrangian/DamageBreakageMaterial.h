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

  virtual void computeQpProperties() override;

  virtual Real updatedamage(); //update damage variable and return 

  virtual void updatebreakage(); //update breakage variable and return

  virtual void updatemodulus(Real alpha_updated); //update modulus

protected:
  
  //declare material properties:
  MaterialProperty<Real> & _alpha_damagedvar; //damage variable
  MaterialProperty<Real> & _B_breakagevar;    //breakage variable

  MaterialProperty<Real> & _shear_modulus;    //shear modulus
  MaterialProperty<Real> & _damaged_modulus;  //damaged modulus
  
  //get old material properties:
  const MaterialProperty<Real> & _alpha_damagedvar_old; //old damage variable
  const MaterialProperty<Real> & _B_breakagevar_old;     //old breakage variable
  const MaterialProperty<Real> & _I2_old;     //old second strain invariants
  const MaterialProperty<Real> & _xi_old;     //old strain invariant ratio
  const MaterialProperty<Real> & _initial_damage; //initial damage 

  //get const values
  Real _shear_modulus_o;
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
  /// coefficient of damage solid modulus
  Real _gamma_damaged_r;
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

};