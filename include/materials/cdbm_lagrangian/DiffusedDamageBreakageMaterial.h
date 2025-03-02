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
class DiffusedDamageBreakageMaterial : public Material
{
public:
  static InputParameters validParams();

  DiffusedDamageBreakageMaterial(const InputParameters & parameters);

  virtual void initQpStatefulProperties() override;

  virtual void computeQpProperties() override;

  //@brief update damage variable and return, Cd here could be strain dependent
  virtual void updateCd(); //update Cd variable and return

  //@brief update C1 C2 variable and return, C1 C2 here could be real or aux variable
  virtual void updateC1C2(); //update C1 C2 variable and return

  //@brief update xi_0 variable and return, xi_0 here could be real or aux variable
  virtual void updatexi0(); //update xi_0 variable and return

  virtual Real updatedamage(); //update damage variable and return 

  virtual void updatebreakage(); //update breakage variable and return

  virtual void updatemodulus(Real alpha_updated); //update modulus

  virtual void computegammar(); //compute gamma_r

  virtual void computecoefficients(); //compute coefficients a_0 a_1 a_2 a_3 of granular phase

  virtual Real alphacr_root1(Real xi); //alpha cr root 1

  virtual Real alphacr_root2(Real xi); //alpha cr root 2

  /**
   * Compute the crack strain in the crack coordinate system. Also
   * computes the crack orientations, and stores in _crack_rotation.
   * @param strain_in_crack_dir Computed strains in crack directions
   */
  void computePrincipalStrainAndOrientation(RealVectorValue & strain_in_crack_dir);  

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

  MaterialProperty<Real> & _xi_1_mat;         //xi

  MaterialProperty<Real> & _xi_0_mat; //xi_0

  //initial shear modulus distribution
  MaterialProperty<Real> & _shear_modulus_o_mat; //shear modulus

  //use "mat" 
  //damage accmulation rate
  MaterialProperty<Real> & _Cd_mat; //Cd_constant

  //damage healing parameters
  MaterialProperty<Real> & _C1_mat; //C1
  MaterialProperty<Real> & _C2_mat; //C2

  //maximum principal strain rate
  MaterialProperty<Real> & _strain_dir0_positive; //strain_dir0_positive_old
  const MaterialProperty<Real> & _strain_dir0_positive_old; //strain_dir0_positive_old

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
  const MaterialProperty<RankTwoTensor> & _elastic_strain_old;   //old first lamé constant
  const MaterialProperty<RankTwoTensor> & _plastic_strain_old;   //old first lamé constant

  //get const values
  Real _lambda_o;
  Real _shear_modulus_o;
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

  /// @brief add option to provide xi0 as a constant value or as an auxiliary variable
  const bool _use_xi0_aux;
  const Real _xi0_value;
  const VariableValue * _xi0_aux;

  /// @brief add option to provide initial shear modulus as a constant value or as an auxiliary variable
  const bool _use_shear_modulus_o_aux;
  const Real _shear_modulus_o_value;
  const VariableValue * _shear_modulus_o_aux;

  /// @brief add option to provide xi as a nonlocal variable
  const bool _use_nonlocal_xi;
  const VariableValue * _nonlocal_xi;  

  /// @brief add option to provide cd as aux variable
  const bool _use_cd_aux;
  const VariableValue * _cd_aux; 

  /// @brief add option to provide cb_multiplier as aux variable
  const bool _use_cb_multiplier_aux;
  const VariableValue * _cb_multiplier_aux;

  /// @brief add option to provide cbh as aux variable
  const bool _use_cbh_aux;
  const VariableValue * _cbh_aux;

  /// @brief add option to provide c1 as aux variable
  const bool _use_c1_aux;
  const VariableValue * _c1_aux;

  /// @brief add option to provide c2 as aux variable
  const bool _use_c2_aux;
  const VariableValue * _c2_aux;

  /// @brief add option to use strain-dependent cd
  const bool _use_cd_strain_dependent;
  // Add member variable for block ID (where the rate-dependent Cd applies)
  unsigned int _block_id;
  const Real _m_exponent;
  const Real _strain_rate_hat;
  const Real _cd_hat;
  const int _block_id_applied; 

  /// @brief add option to use plastic strain rate
  // default is to use elastic strain rate
  const bool _use_plastic_strain_rate;

  /// @brief add damage variable value from system solve
  const VariableValue & _alpha_damagedvar_sys;
};