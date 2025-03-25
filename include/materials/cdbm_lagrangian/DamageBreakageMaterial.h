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

  virtual void computecoefficientsgivenxi(); //compute coefficients a_0 a_1 a_2 a_3 with given xi for energy difference

  // virtual void computecoefficientsgivenchi(); //compute coefficients a_0 a_1 a_2 a_3 with given chi for energy difference

  //2D only
  virtual void computeinitialdamage2D(); //compute initial damage with time dependent material

  //3D only
  virtual void computeinitialdamage3D(); //compute initial damage with time dependent material

  //2D only
  virtual void computedamageperturbation2D(); //compute damage perturbation with time dependent material

  //3D only
  virtual void computedamageperturbation3D(); //compute damage perturbation with time dependent material

  virtual void buildLmatrix(); //build L matrix

  virtual Real alphacr_root1(Real xi); //alpha cr root 1

  virtual Real alphacr_root2(Real xi); //alpha cr root 2

  /* Additional Features */
  virtual void addDilatancyCompactionAnand(); //add dilatancy/compaction effect using anand model

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

  MaterialProperty<Real> & _initial_damage_time_dependent_mat; //initial damage handling time dependent damage
  MaterialProperty<Real> & _damage_perturbation; //damage perturbation
  const MaterialProperty<Real> & _damage_perturbation_old; //old damage perturbation

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

  //damage accmulation rate
  MaterialProperty<Real> & _Cd_rate_dependent; //Cd_constant

  //maximum principal strain rate
  MaterialProperty<Real> & _strain_dir0_positive; //strain_dir0_positive_old
  const MaterialProperty<Real> & _strain_dir0_positive_old; //strain_dir0_positive_old

  //Velocity gradient L
  MaterialProperty<bool> & _use_vels_build_L_mat; //use_vels_build_L
  MaterialProperty<RankTwoTensor> & _velgrad_L; //L

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
  const MaterialProperty<RankTwoTensor> & _elastic_strain_old;   //old elastic strain
  const MaterialProperty<RankTwoTensor> & _total_lagrange_strain_old;   //old total lagrangian strain

  //use velocity to build L
  const bool _use_vels_build_L;
  const VariableGradient * _grad_vel_x;
  const VariableGradient * _grad_vel_y;
  const VariableGradient * _grad_vel_z;

  //get const values
  Real _lambda_o;
  Real _shear_modulus_o;
  /// strain invariants ratio: onset of breakage healing
  Real _xi_d;
  /// strain invariants ratio: minimum allowable value
  Real _xi_min;
  /// strain invariants ratio: maximum allowable value
  Real _xi_max;
  /// strain invariants ratio: given value
  Real _xi_given;
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

  /// @brief add option to provide cg as aux variable
  const bool _use_cg_aux;
  const VariableValue * _cg_aux;

  /// @brief add option to provide xi as aux variable
  const bool _use_const_xi_aux;
  const VariableValue * _const_xi_aux;
  const int _const_xi_block_id;

  /// @brief add option to use strain-dependent cd
  const bool _use_cd_strain_dependent;
  // Add member variable for block ID (where the rate-dependent Cd applies)
  unsigned int _block_id;
  const Real _m_exponent;
  const Real _strain_rate_hat;
  const Real _cd_hat;
  const int _straindep_block_id_applied; 

  /// @brief add option to use total strain rate
  // default is to use elastic strain rate
  const bool _use_total_strain_rate;

  /// @brief add option to use pore pressure to decrease mean stress
  // pore pressure is saved as material property and will be used 
  // in ComputeLagrangianDamageBreakageStressPK2.C
  const bool _use_pore_pressure;
  const VariableValue * _pore_pressure; 
  MaterialProperty<Real> & _pore_pressure_mat;

  /// @brief add option to use overstress to nucleate the rupture
  const bool _use_overstress;
  const VariableValue * _overstress;
  MaterialProperty<bool> & _use_overstress_mat;
  MaterialProperty<Real> & _overstress_mat;

  /// @brief add option to build initial damage profile inside this material object
  const bool _build_param_use_initial_damage_time_dependent_mat;
  const Real _build_param_peak_value;
  const Real _build_param_sigma;
  const Real _build_param_len_of_fault;

  /// @brief add option to build initial damage profile in 3D
  const bool _build_param_use_initial_damage_3D;
  const Real _build_param_len_of_fault_dip;
  const std::vector<Real> _build_param_center_point;

  /// @brief add option to build time dependent damage perturbation inside this material object
  const bool _perturbation_build_param_use_damage_perturb;
  const std::vector<Real> _perturbation_build_param_nucl_center;
  const Real _perturbation_build_param_length;
  const Real _perturbation_build_param_peak_value;
  const Real _perturbation_build_param_sigma;
  const Real _perturbation_build_param_thickness;
  const Real _perturbation_build_param_duration;

  /// @brief add constant parameters for state variable evolution
  const bool _use_state_var_evolution;
  const Real _const_A;
  const Real _const_B;
  const Real _const_theta_o;
  const Real _initial_theta0;
  MaterialProperty<bool> & _use_state_var_evolution_mat;
  MaterialProperty<Real> & _const_A_mat;
  MaterialProperty<Real> & _const_B_mat;
  MaterialProperty<Real> & _const_theta_o_mat;
  MaterialProperty<Real> & _initial_theta0_mat;

  /// @brief add option to use energy breakage evolution equation
  const bool _use_energy_breakage_evolution;

  /// @brief add option to add dilatancy/compaction effect
  const bool _add_dilatancy_compaction_anand;
  const Real _anand_param_go;
  const Real _anand_param_eta_cv;
  const Real _anand_param_p;
  MaterialProperty<bool> & _add_dilatancy_compaction_anand_mat;
  MaterialProperty<Real> & _anand_param_go_mat;
  MaterialProperty<Real> & _anand_param_eta_cv_mat;
  MaterialProperty<Real> & _anand_param_p_mat;

};