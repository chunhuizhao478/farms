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
 * ComputePoroDamageBreakageStress3DSlipWeakeningnonlocal put everything inside the computeQpstress without defining
 * additional functions
 
 */
class ComputePoroDamageBreakageStress3DSlipWeakeningnonlocal : public ComputeDamageBreakageStressBase3D
{
public:
  static InputParameters validParams();

  ComputePoroDamageBreakageStress3DSlipWeakeningnonlocal(const InputParameters & parameters);

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
  virtual void computeQpTangentModulus(RankFourTensor & tangent, Real I1, Real I2, Real xi, RankTwoTensor Ee, 
                                       Real a0, Real a1, Real a2, Real a3, Real gamma_damaged_r);

  /// @brief Setup initial values for the first step
  void setupInitial();

  /// @brief Compute deviatoric stress tensor
  void computeDeviatroicStrainRateTensor();

  /// @brief Compute strain rate dependent Cd
  void computeStrainRateCd();

  void computeSimplifiedPlasticDerivatives(const RankTwoTensor & eps_p,
                                         const RankTwoTensor & eps_p_dot);

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
  const MaterialProperty<RankTwoTensor> & _sts_total_old;

  /// Get initial values
  const MaterialProperty<RankTwoTensor> & _static_initial_stress_tensor;
  const MaterialProperty<RankTwoTensor> & _static_initial_strain_tensor;
  const MaterialProperty<RankTwoTensor> & _sts_initial_tensor_old;

  const MaterialProperty<Real> & _initial_porepressure;

  // const MaterialProperty<Real> & _I1_initial;
  // const MaterialProperty<Real> & _I2_initial;
  // const MaterialProperty<Real> & _xi_initial;
  const MaterialProperty<Real> & _initial_damage;
  const MaterialProperty<Real> & _initial_breakage;

  /// perturbation (damage)
  const MaterialProperty<Real> & _damage_perturbation;
  /// perturbation (shear stress)
  // const MaterialProperty<Real> & _shear_stress_perturbation;

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

  int & _step;

  /// matprop : deviatoric stress tensor
  MaterialProperty<Real> & _deviatroic_strain_rate;
  const MaterialProperty<Real> & _deviatroic_strain_rate_old;

  /// matprop : Cd
  MaterialProperty<Real> & _Cd_mat;
  const MaterialProperty<Real> & _Cd_mat_old;

  /// strain rate dependent Cd parameters
  bool _use_strain_rate_dependent_Cd;
  Real _m_exponent;
  Real _strain_rate_hat;
  Real _cd_hat;

  /// option: set Cd = 0 when strain rate < strain_rate_hat (default: false -> use cd_hat)
  bool _zero_Cd_below_threshold;

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

  /// nonlocal equivalent strain
  bool _use_nonlocal_eqstrain;
  const MaterialProperty<Real> & _eqstrain_nonlocal_old;

  /// blocks where nonlocal equivalent strain is enabled; empty means all blocks
  const std::vector<unsigned int> _nonlocal_eqstrain_blocks;

  /// nonlocal strain rate for Cd calculation
  bool _use_nonlocal_strain_rate;
  const MaterialProperty<Real> * _strain_rate_nonlocal_old;

  /// static solve flag
  bool _static_solve_flag;

  /// helper: whether nonlocal eqstrain should be used on the current element
  inline bool useNonlocalEqStrainHere() const
  {
    if (!_use_nonlocal_eqstrain)
      return false;
    if (_nonlocal_eqstrain_blocks.empty())
      return true;
  const auto sid = static_cast<unsigned int>(_current_elem->subdomain_id());
    return std::find(_nonlocal_eqstrain_blocks.begin(), _nonlocal_eqstrain_blocks.end(), sid) !=
           _nonlocal_eqstrain_blocks.end();
  }
  
};