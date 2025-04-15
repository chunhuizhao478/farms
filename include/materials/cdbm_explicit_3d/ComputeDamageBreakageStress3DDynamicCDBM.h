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
 * ComputeDamageBreakageStress3DDynamicCDBM put everything inside the computeQpstress without defining
 * additional functions
 
 */
class ComputeDamageBreakageStress3DDynamicCDBM : public ComputeDamageBreakageStressBase3D
{
public:
  static InputParameters validParams();

  ComputeDamageBreakageStress3DDynamicCDBM(const InputParameters & parameters);

  virtual void initialSetup() override;

protected:
  virtual void initQpStatefulProperties() override;
  virtual void computeQpStress() override;

  /// @brief Compute gamma_r
  /// @return gamma_r
  Real computegammar();

  /// @brief Compute breakage coefficients
  /// @param gamma_damaged_r
  /// @return a0 a1 a2 a3 in a vector
  std::vector<Real> computecoefficients(Real gamma_damaged_r);

  /// @brief Compute first root of hessian matrix
  /// @param xi 
  /// @return the first root of critical alpha_cr
  Real alphacr_root1(Real xi, Real gamma_damaged_r);

  /// @brief Compute second root of hessian matrix
  /// @param xi 
  /// @return the second root of critical alpha_cr
  Real alphacr_root2(Real xi, Real gamma_damaged_r);

  /// @brief Compute elasticity tensor for small strain
  virtual void computeQpTangentModulus(RankFourTensor & tangent, Real I1, Real I2, Real xi, RankTwoTensor Ee);

  /**
   * Compute the crack strain in the crack coordinate system. Also
   * computes the crack orientations, and stores in _crack_rotation.
   * @param strain_in_crack_dir Computed strains in crack directions
   */
  void computeCrackStrainAndOrientation(RealVectorValue & strain_in_crack_dir);
  
  /**
   * Update the solid bulk compliance due to the effects of cracking
   */
  void updateSolidBulkCompliance();

  /**
   *  Update the effective permeability tensor to account for the effect of cracking
   */
  void updatePermeabilityForCracking();  

  /**
   *  Implement the Heaviside function
   */
  void HeavisideFunction(Real & x);  

  /**
   *  Implement the effective stress modification
   */
  void computeEffectiveStressModification();

  /// additional variables
  /// strain invariants ratio: onset of damage evolution
  Real _xi_0;

  /// strain invariants ratio: onset of breakage healing
  Real _xi_d;

  /// strain invariants ratio: minimum allowable value
  Real _xi_min;

  /// strain invariants ratio: maximum allowable value
  Real _xi_max;

  /// energy ratio
  Real _chi;

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
  const MaterialProperty<RankTwoTensor> & _crack_rotation_old;

  //add grad term
  const VariableValue & _alpha_grad_x;
  const VariableValue & _alpha_grad_y;
  const VariableValue & _alpha_grad_z;

  /// diffusion coefficient
  Real _D;

  /// Get initial values
  // const MaterialProperty<RankTwoTensor> & _static_initial_stress_tensor;
  // const MaterialProperty<RankTwoTensor> & _static_initial_strain_tensor;
  // const MaterialProperty<Real> & _I1_initial;
  // const MaterialProperty<Real> & _I2_initial;
  // const MaterialProperty<Real> & _xi_initial;
  const MaterialProperty<Real> & _initial_damage;
  const MaterialProperty<Real> & _initial_breakage;

  /// shear stress perturbation
  const MaterialProperty<Real> & _shear_stress_perturbation;

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

  /// dimension
  const unsigned int _dim;

  /// @brief include porous flow coupling
  /// @return true if porous flow coupling is included
  bool _porous_flow_coupling;

  /// @brief define the crack surface roughness correction factor
  const Real _crack_surface_roughness_correction_factor;

  /// @brief define the length scale
  const Real _length_scale;

  /// @brief define the intrisic permeability
  const Real _intrinsic_permeability; 

  /// @brief effective stress modification
  /// @return true if effective stress modification is included
  bool _effective_stress_modification;
  /// @brief fluid density
  Real _fluid_density;
  /// @brief gravity
  Real _gravity;

  //## DEBUG ##
  /// Name of the elasticity tensor material property
  const std::string _elasticity_tensor_name;
  /// Elasticity tensor material property
  const MaterialProperty<RankFourTensor> & _elasticity_tensor;

};