//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ComputeStressBase.h"

/**
 * ComputeDamageBreakageStress3DDynamicCDBMDiffused put everything inside the computeQpstress without defining
 * additional functions
 
 */
class ComputeDamageBreakageStress3DDynamicCDBMDiffused : public ComputeStressBase
{
public:
  static InputParameters validParams();

  ComputeDamageBreakageStress3DDynamicCDBMDiffused(const InputParameters & parameters);

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
  virtual void computeQpTangentModulus(RankFourTensor & tangent, Real I1, Real I2, Real xi, RankTwoTensor Ee, Real a0, Real a1, Real a2, Real a3, Real gamma_damaged_r);

  /// @brief Compute deviatoric stress tensor
  void computeDeviatroicStrainRateTensor();

  /// @brief Compute strain rate dependent Cd
  void computeStrainRateCd();

  void computeSimplifiedPlasticDerivatives(const RankTwoTensor & eps_p,
                                        const RankTwoTensor & eps_p_dot);

  /// additional variables
  /// initial lambda constant value
  Real _lambda_o;

  /// initial shear modulus value
  Real _shear_modulus_o;

  /// strain invariants ratio: onset of damage evolution
  Real _xi_0;

  /// strain invariants ratio: onset of breakage healing
  Real _xi_d;

  /// energy ratio
  Real _chi;

  /// material parameter: compliance or fluidity of the fine grain granular material
  Real _C_g;

  /// coefficient of power law indexes
  Real _m1;

  /// coefficient of power law indexes
  Real _m2;

  /// alpha_damagedvar_aux
  const VariableValue & _alpha_damagedvar_aux;

  /// B_damagedvar_aux
  const VariableValue & _B_damagedvar_aux;

  /// deviatoric strain rate
  MaterialProperty<Real> & _deviatroic_strain_rate;

  /// deviatroic stress
  MaterialProperty<RankTwoTensor> & _sigma_d;
  const MaterialProperty<RankTwoTensor> & _sigma_d_old;

  /// plastic strain
  MaterialProperty<RankTwoTensor> & _eps_p;
  const MaterialProperty<RankTwoTensor> & _eps_p_old;

  /// total strain
  MaterialProperty<RankTwoTensor> & _eps_total;
  const MaterialProperty<RankTwoTensor> & _eps_total_old;
  
  /// elastic strain
  MaterialProperty<RankTwoTensor> & _eps_e;

  /// I1
  MaterialProperty<Real> & _I1;

  /// I2
  MaterialProperty<Real> & _I2;

  /// xi
  MaterialProperty<Real> & _xi;

  /// lambda
  MaterialProperty<Real> & _lambda;

  /// shear_modulus
  MaterialProperty<Real> & _shear_modulus;

  /// damaged_modulus
  MaterialProperty<Real> & _gamma_damaged;

  /// shear stress perturbation
  const MaterialProperty<Real> & _stress_perturbation;

   /// Biot coefficient for solid phase
  MaterialProperty<Real> & _Biot_coeff_s;

  /// Biot coefficient for granular phase
  MaterialProperty<Real> & _Biot_coeff_g;

  /// Biot modulus for solid phase
  MaterialProperty<Real> & _Biot_modulus_s;

  /// Biot modulus for granular phase
  MaterialProperty<Real> & _Biot_modulus_g;

  /// effective Biot coefficient
  MaterialProperty<Real> & _biot_coeff_eff;

  /// effective Biot modulus
  MaterialProperty<Real> & _Biot_modulus_eff;

  /// fluid-solid coupling parameter
  MaterialProperty<Real> & _fluid_solid_coupling;

  /// permeability of solid phase
  MaterialProperty<Real> & _perm_s;

  /// permeability of granular phase
  MaterialProperty<Real> & _perm_g;

  /// critical permeability
  MaterialProperty<Real> & _perm_cr;

  /// critical porosity
  MaterialProperty<Real> & _phi_cr;

  /// plastic porosity
  MaterialProperty<Real> & _phi_p;

  /// initial permeability of solid material
  Real _permeability_solid_o;

  /// solid bulk modulus of solid grains
  Real _solid_bulk_modulus_s;

  /// solid bulk modulus of granular material
  Real _solid_bulk_modulus_g;

  /// fluid bulk modulus
  Real _fluid_bulk_modulus;

  /// initial porosity of solid phase
  Real _porosity_solid_o;

  /// initial fluid viscosity
  Real _initial_viscosity_fluid;

  /// parameter for permeability evolution with damage
  Real _b;

  /// initial harmonic mean grain size
  Real _DHo;

  /// ultimate harmonic mean grain size
  Real _DHu;

  /// pore pressure coupled variable
  const VariableValue & _pore_pressure;

  /// old pore pressure coupled variable
  const VariableValue & _pore_pressure_old;

  /// PK1 off-diagonal Jacobian
  MaterialProperty<RankTwoTensor> & _stress_off_diag_jacobian;

  /// use dilatancy
  bool _use_dilatancy;

  /// plastic volume change
  MaterialProperty<Real> & _eta;

  /// old plastic volume change
  const MaterialProperty<Real> & _eta_old;

  /// dilatancy function beta
  MaterialProperty<Real> & _dilatancy_function_beta;

  /// shear rate nu tensor
  MaterialProperty<RankTwoTensor> & _shear_rate_nu;

  /// Anand parameter go
  Real _anand_param_go_mat;

  /// Anand parameter eta_cv
  Real _anand_param_eta_cv_mat;

  /// Anand parameter p
  Real _anand_param_p_mat;

  MaterialProperty<RankTwoTensor> & _deps_p_dp;      // ∂εᵖ/∂p (simplified FD)
  MaterialProperty<RankFourTensor> & _deps_p_deps;   // ∂εᵖ/∂ε (simplified FD)

  const MaterialProperty<RankTwoTensor> & _mechanical_strain_old;
};