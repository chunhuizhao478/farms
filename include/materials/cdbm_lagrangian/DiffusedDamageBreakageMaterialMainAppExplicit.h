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
class DiffusedDamageBreakageMaterialMainAppExplicit : public Material
{
public:
  static InputParameters validParams();

  DiffusedDamageBreakageMaterialMainAppExplicit(const InputParameters & parameters);

  virtual void initQpStatefulProperties() override;

  virtual void computeQpProperties() override;

  virtual void computegammar(); //compute gamma_r
  virtual void updatedamagebreakage(); //update damage variable and breakage variable
  virtual void updatemodulus(); //update modulus
  virtual void updateporosolid(); //update poroelastic proeperties for solid material
  virtual void updateporogranular(); //update poroelastic proeperties for granular material
  virtual void computecoefficients(); //compute coefficients: a0 a1 a2 a3
  virtual Real alphacr_root1(Real xi); //compute alpha_cr
  virtual void buildLmatrix(); //build L matrix
  virtual void acceptspatialCg(); //accept spatial cg
  virtual void usestatevar(); //use state variable evolution

protected:
  
  //declare material properties:
  MaterialProperty<Real> & _gamma_damaged_r; //gamma_r
  MaterialProperty<Real> & _alpha_damagedvar; //alpha_damagedvar
  MaterialProperty<Real> & _B_damagedvar; //B_damagedvar
  MaterialProperty<Real> & _lambda; //lambda
  MaterialProperty<Real> & _shear_modulus; //shear_modulus
  MaterialProperty<Real> & _damaged_modulus; //damaged_modulus
  MaterialProperty<Real> & _a0; //a0
  MaterialProperty<Real> & _a1; //a1
  MaterialProperty<Real> & _a2; //a2
  MaterialProperty<Real> & _a3; //a3
  MaterialProperty<Real> & _Biot_coeff_s; //Biot_coefficient_solid
  MaterialProperty<Real> & _Biot_coeff_g; //Biot_coefficient_granular
  MaterialProperty<Real> & _Biot_modulus_s; //Biot_modulus_solid
  MaterialProperty<Real> & _Biot_modulus_g; //Biot_modulus_granular
  MaterialProperty<Real> & _biot_coeff_eff; //_biot_coeff_eff
  MaterialProperty<Real> & _Biot_modulus_eff; //Biot_modulus_eff
  MaterialProperty<Real> & _fluid_solid_coupling; //Fluid_solid_coupling term
  MaterialProperty<Real> & _perm_s; //permeability_solid
  MaterialProperty<Real> & _perm_g; //permeability_granular
  MaterialProperty<Real> & _perm_cr; //permeability_critical
  MaterialProperty<Real> & _phi_cr; //critical porosity
  MaterialProperty<Real> & _phi_p; // plastic porosity
  MaterialProperty<Real> & _fluid_viscosity; // fluid_viscosity
  MaterialProperty<Real> & _C_g; //Cg
  MaterialProperty<Real> & _m1; //m1
  MaterialProperty<Real> & _m2; //m2
  MaterialProperty<Real> & _structural_stress_coefficient; //structral_stress_coefficient
  MaterialProperty<RealGradient> & _grad_alpha_damagedvar; //gradient_alpha_damagedvar
  MaterialProperty<Real> & _grad_alpha_damagedvar_xdir; //gradient_alpha_damagedvar_xdir
  MaterialProperty<Real> & _grad_alpha_damagedvar_ydir; //gradient_alpha_damagedvar_ydir
  MaterialProperty<RankTwoTensor> & _velgrad_L; //L

  //input values
  Real _lambda_o_value; //lambda_o
  Real _shear_modulus_o_value; //shear_modulus_o
  Real _xi_0_value; //xi_0
  Real _xi_d_value; //xi_d
  Real _chi_value; //chi 
  Real _C_g_value; //Cg
  Real _m1_value; //m1
  Real _m2_value; //m2
  Real _permeability_solid_o; //permeability_solid_o
  Real _solid_bulk_modulus_s; //solid_bulk_modulus_solid_o
  Real _solid_bulk_modulus_g; //solid_bulk_modulus_granular_o
  Real _fluid_bulk_modulus; //fluid_bulk_modulus
  Real _porosity_solid_o; //initial prosoity of solid phase
  Real _initial_viscosity_fluid; //initial _initial_viscosity_fluid
  Real _b; //permeability evolution with damage
  Real _DHo; //initial grain size
  Real _DHu; //ultimate grain size

  //input coupled variables from main app
  const VariableValue & _alpha_damagedvar_aux; //alpha_damagedvar
  const VariableValue & _B_damagedvar_aux; //B_damagedvar

  //get structural stress coefficient and gradient of damage variable
  const VariableValue & _structural_stress_coefficient_aux; //structral_stress_coefficient 
  const VariableGradient & _grad_alpha_damagedvar_value; //gradient_alpha_damagedvar

  //use velocity to build L matrix
  const VariableGradient & _grad_vel_x;
  const VariableGradient & _grad_vel_y;
  const VariableGradient & _grad_vel_z;

  /// First Elastic Strain Invariant
  const MaterialProperty<Real> & _I1;
  /// Elastic deformation gradient
  const MaterialProperty<RankTwoTensor> & _Fe;
  /// Second Elastic Strain Invariant
  const MaterialProperty<Real> & _xi;
  /// Plastic Jacobian
  const MaterialProperty<Real> & _Jp;
  /// Plastic strain rate tensor
  const MaterialProperty<RankTwoTensor> & _Dp_old;
  /// Plastic strain rate older tensor
  const MaterialProperty<RankTwoTensor> & _Dp_older;
    /// Plastic porosity old
  const MaterialProperty<Real> & _phi_p_old;
  /// Plastic porosity older
  const MaterialProperty<Real> & _phi_p_older;

  //use spatial cg
  bool _use_spatial_cg; //use spatial cg
  const VariableValue & _cg_aux; //cg_aux

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

};