//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ComputeLagrangianDamageBreakageStressPK2Diffused.h"


registerMooseObject("farmsApp", ComputeLagrangianDamageBreakageStressPK2Diffused);

InputParameters
ComputeLagrangianDamageBreakageStressPK2Diffused::validParams()
{
  InputParameters params = ComputeLagrangianStressPK1::validParams();
  // Add pore pressure gradient coupling
  params.addRequiredCoupledVar(
    "porepressure", 
    "The pore pressure appropriate for the simulation geometry and coordinate system");
  
  //use dilatancy dependent variables
  params.addParam<bool>("use_dilatancy", false, "Flag to use dilatancy variable evolution");
  params.addParam<Real>("anand_param_go_mat",0.2,"Dilatancy parameter go");
  params.addParam<Real>("anand_param_eta_cv_mat",0.006,"Dilatancy parameter eta_cv");
  params.addParam<Real>("anand_param_p_mat",2,"Dilatancy parameter p");

  return params;
}

ComputeLagrangianDamageBreakageStressPK2Diffused::ComputeLagrangianDamageBreakageStressPK2Diffused(const InputParameters & parameters)
  : ComputeLagrangianStressPK1(parameters),
  _pore_pressure(coupledValue("porepressure")),
  _pore_pressure_old(coupledValueOld("porepressure")),
  _pk1_off_diag_jacobian(declareProperty<RankTwoTensor>(_base_name + "pk1_off_diag_jacobian")),
  _dI1dF(declareProperty<RankTwoTensor>(_base_name + "first_elastic_strain_invariant_derivative")),
  _dJpdF(declareProperty<RankTwoTensor>(_base_name + "plastic_jacobian_derivative")),
  _dJpdp(declareProperty<Real>(_base_name + "plastic_jacobian_derivative_pressure")),
  _dDpdF(declareProperty<RankFourTensor>(_base_name + "plastic_deformation_rate_derivative")),
  _dDpdp(declareProperty<RankTwoTensor>(_base_name + "plastic_deformation_rate_derivative_pressure")),
  _Fp(declareProperty<RankTwoTensor>(_base_name + "plastic_deformation_gradient")),
  _Jp(declareProperty<Real>(_base_name + "plastic_deformation_gradient_det")),
  _Fe(declareProperty<RankTwoTensor>(_base_name + "elastic_deformation_gradient")),
  _Tau(declareProperty<RankTwoTensor>(_base_name + "deviatroic_stress")),
  _Ee(declareProperty<RankTwoTensor>(_base_name + "green_lagrange_elastic_strain")),
  _Ep(declareProperty<RankTwoTensor>(_base_name + "plastic_strain")),
  _E(declareProperty<RankTwoTensor>(_base_name + "total_lagrange_strain")),
  _I1(declareProperty<Real>(_base_name + "first_elastic_strain_invariant")),
  _I2(declareProperty<Real>(_base_name + "second_elastic_strain_invariant")),
  _xi(declareProperty<Real>(_base_name + "strain_invariant_ratio")),
  _S(declareProperty<RankTwoTensor>(_base_name + "pk2_stress")),  
  _Tp(declareProperty<RankTwoTensor>(_base_name + "plastic_stress")),
  _C(declareProperty<RankFourTensor>(_base_name + "pk2_jacobian")),
  _Dp(declareProperty<RankTwoTensor>(_base_name + "plastic_strain_rate")),
  _Fp_dot(declareProperty<RankTwoTensor>(_base_name + "cdbm_plastic_deformation_gradient_rate")),
  _F_dot(declareProperty<RankTwoTensor>(_base_name + "cdbm_deformation_gradient_rate")),
  _D(declareProperty<RankTwoTensor>(_base_name + "deformation_rate")),
  //---------------------------------------------------------------------------------------------//
  // Add option to add dilatancy/compaction effect //Follow paper Section 7.1
  _use_dilatancy(getParam<bool>("use_dilatancy")),
  _eta(declareProperty<Real>(_base_name + "plastic_volume_change")),
  _eta_old(getMaterialPropertyOldByName<Real>(_base_name + "plastic_volume_change")),
  _dilatancy_function_beta(declareProperty<Real>(_base_name + "dilatancy_function_beta")),
  _shear_rate_nu(declareProperty<RankTwoTensor>(_base_name + "shear_rate_nu")),
  _anand_param_go_mat(getParam<Real>("anand_param_go_mat")),
  _anand_param_eta_cv_mat(getParam<Real>("anand_param_eta_cv_mat")),
  _anand_param_p_mat(getParam<Real>("anand_param_p_mat")),
  //---------------------------------------------------------------------------------------------//
  _deviatroic_strain_rate(declareProperty<Real>("deviatroic_strain_rate")),
  //---------------------------------------------------------------------------------------------//
  _lambda_const(getMaterialProperty<Real>("lambda_const")),
  _shear_modulus(getMaterialProperty<Real>("shear_modulus")),
  _damaged_modulus(getMaterialProperty<Real>("damaged_modulus")),
  _B_breakagevar(getMaterialProperty<Real>("B_damagedvar")),
  _B_breakagevar_old(getMaterialPropertyOldByName<Real>("B_damagedvar")),
  _Tau_old(getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "deviatroic_stress")),
  _Fp_old(getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "plastic_deformation_gradient")),
  _F_old(getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "deformation_gradient")),
  _Ep_old(getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "plastic_strain")),
  //---------------------------------------------------------------------------------------------//
  _C_g(getMaterialProperty<Real>("C_g")),
  _m1(getMaterialProperty<Real>("m1")),
  _m2(getMaterialProperty<Real>("m2")),
  _a0(getMaterialProperty<Real>("a0")),
  _a1(getMaterialProperty<Real>("a1")),
  _a2(getMaterialProperty<Real>("a2")),
  _a3(getMaterialProperty<Real>("a3")),
  //---------------------------------------------------------------------------------------------//
  // Poroelastic properties
  _Biot_coeff_s(getMaterialProperty<Real>("Biot_coefficient_solid")),
  _Biot_coeff_g(getMaterialProperty<Real>("Biot_coefficient_granular")),
  _Biot_modulus_s(getMaterialProperty<Real>("Biot_modulus_solid")),
  _Biot_modulus_g(getMaterialProperty<Real>("Biot_modulus_granular")),
  //---------------------------------------------------------------------------------------------//
  _structural_stress_coefficient(getMaterialProperty<Real>("structural_stress_coefficient")),
  _grad_alpha_damagedvar(getMaterialProperty<RealGradient>("gradient_alpha_damagedvar")),
  //---------------------------------------------------------------------------------------------//
  //For the strain rate Cd
  _velgrad_L(getMaterialProperty<RankTwoTensor>("velgrad_L")),
  //add shear stress perturbation
  //_dim(_mesh.dimension())
  //---------------------------------------------------------------------------------------------//
  //use state dependent variables
  _Theta(declareProperty<Real>("state_variable")),
  _Theta_old(getMaterialPropertyOldByName<Real>("state_variable")),
  _use_state_var_evolution_mat(getMaterialProperty<bool>("use_state_var_evolution_mat")),
  _const_A_mat(getMaterialProperty<Real>("const_A_mat")),
  _const_B_mat(getMaterialProperty<Real>("const_B_mat")),
  _const_theta_o_mat(getMaterialProperty<Real>("const_theta_o_mat")),
  _initial_theta0_mat(getMaterialProperty<Real>("initial_theta0_mat")),
  //---------------------------------------------------------------------------------------------//
  //add shear stress perturbation
  _shear_stress_perturbation(getMaterialPropertyOldByName<Real>("shear_stress_perturbation"))
{
}

//Rules:See https://github.com/idaholab/moose/discussions/19450
//Only the object that declares the material property can assign values to it.
//Objects can request material properties, gaining read-only access to their values.
//When any object (including the object that declares it) requests the old value of a material property, that property becomes "stateful".
//All stateful material properties must be initialized within the initQpStatefulProperties call. 
//
void
ComputeLagrangianDamageBreakageStressPK2Diffused::initQpStatefulProperties()
{
  _Fp[_qp] = RankTwoTensor::Identity();
  _Fe[_qp] = RankTwoTensor::Identity();
  _Tau[_qp].zero();
  _Ee[_qp].zero();
  _Ep[_qp].zero();
  _E[_qp].zero();
  _I1[_qp] = 0.0;
  _I2[_qp] = 0.0;
  _xi[_qp] = -sqrt(3);
  _S[_qp].zero();
  _C[_qp].zero();
  _eta[_qp] = 0.0;
  _dilatancy_function_beta[_qp] = 0.0;
  _shear_rate_nu[_qp] = 0.0;
  _Theta[_qp] = _initial_theta0_mat[_qp];
  // Initialize PK1 off-diagonal Jacobian
  _pk1_off_diag_jacobian[_qp].zero();
  _dI1dF[_qp].zero();
  _dJpdF[_qp].zero();
  _dJpdp[_qp] = 0.0;
  _dDpdF[_qp].zero();
  _dDpdp[_qp].zero();

}

void
ComputeLagrangianDamageBreakageStressPK2Diffused::computeQpPK1Stress()
{

  //--------------------------------------------------------------------------
  // PK2 update
  computeQpPK2Stress();

  _Jp[_qp] = _Fp[_qp].det();
  
  RankTwoTensor Fpinv = _Fp[_qp].inverse();

  // Compute Fp_dot, F_dot
  // Here we approximate the rate by first-order, not sure if this is sufficient for varying time steps
  // currently MOOSE don't support getMaterialPropertyDot
  RankTwoTensor Fp_dot; Fp_dot.zero();
  RankTwoTensor F_dot; F_dot.zero();

  for (unsigned int i = 0; i < 3; i++){
    for (unsigned int j = 0; j < 3; j++){
        //F_dot(i,j)  = (_F[_qp](i,j) - _F_old[_qp](i,j) ) / _dt; 
        F_dot(i,j)  = (_F[_qp](i,j) - _F_old[_qp](i,j) ); 
        for (unsigned int m = 0; m < 3; m++){
          Fp_dot(i,j) += _Dp[_qp](i,m) * _Fp[_qp](m,j);
        }
    }
  }

  //--------------------------------------------------------------------------
  // Precompute the 4D tensor dFpdF_tensor[i][j][k][l] = dFpdF(i,j,k,l)
  // where dFpdF(i,j,k,l) = 0 if _dt==0 or if either Fp_dot(i,j) or F_dot(k,l) vanish;
  // otherwise dFpdF = Fp_dot(i,j)/F_dot(k,l)
  Real dFpdF_tensor[3][3][3][3];
  for (unsigned int i = 0; i < 3; i++){
    for (unsigned int j = 0; j < 3; j++){
      for (unsigned int k = 0; k < 3; k++){
        for (unsigned int l = 0; l < 3; l++){
          if (_dt == 0.0 || _Fp_dot[_qp](i,j) == 0.0 || _F_dot[_qp](k,l) == 0.0)
            dFpdF_tensor[i][j][k][l] = 0.0;
          else
            dFpdF_tensor[i][j][k][l] = _Fp_dot[_qp](i,j) / _F_dot[_qp](k,l);
        }
      }
    }
  }

  // Define a simple inline delta function
  auto delta = [](int i, int j) -> Real { return (i == j) ? 1.0 : 0.0; };

  //--------------------------------------------------------------------------
  // Precompute dJpdF(k,l)
  // dJpdF(k,l) = Jp * Fpinv(n,m) * dFpdF(m,n,k,l)
  Real dJpdF_tensor[3][3];
  for (unsigned int k = 0; k < 3; k++){
    for (unsigned int l = 0; l < 3; l++){
      Real sum = 0.0;
      for (unsigned int m = 0; m < 3; m++){
        for (unsigned int n = 0; n < 3; n++){
          sum += _Jp[_qp] * Fpinv(n, m) * dFpdF_tensor[m][n][k][l];
        }
      }
      dJpdF_tensor[k][l] = sum;
    }
  }

  //--------------------------------------------------------------------------
  // Precompute dFedF_tensor(i, m, k, l)
  // dFedF(i,m,k,l) = delta(i,k)*Fpinv(l,m) - sum_{h,r}[ _Fe(i,h)*dFpdF(h,r,k,l)*Fpinv(r,m) ]
  Real dFedF_tensor[3][3][3][3];
  for (unsigned int i = 0; i < 3; i++){
    for (unsigned int m = 0; m < 3; m++){
      for (unsigned int k = 0; k < 3; k++){
        for (unsigned int l = 0; l < 3; l++){
          Real val = delta(i, k) * Fpinv(l, m);
          for (unsigned int h = 0; h < 3; h++){ //summation applies to h r
            for (unsigned int r = 0; r < 3; r++){
              val -= _Fe[_qp](i, h) * dFpdF_tensor[h][r][k][l] * Fpinv(r, m);
            }
          }
          dFedF_tensor[i][m][k][l] = val;
        }
      }
    }
  }

  //--------------------------------------------------------------------------
  // Precompute dEdF_tensor(p,q,k,l)
  // dEdF(p,q,k,l) = 0.5 * sum_{m}[ dFedF(m,p,k,l)*_Fe(m,q) + _Fe(m,p)*dFedF(m,q,k,l) ]
  Real dEdF_tensor[3][3][3][3];
  for (unsigned int p = 0; p < 3; p++){
    for (unsigned int q = 0; q < 3; q++){
      for (unsigned int k = 0; k < 3; k++){
        for (unsigned int l = 0; l < 3; l++){
          Real sum = 0.0;
          for (unsigned int m = 0; m < 3; m++){ //summation applies to m
            sum += 0.5 * ( dFedF_tensor[m][p][k][l] * _Fe[_qp](m, q)
                         + _Fe[_qp](m, p) * dFedF_tensor[m][q][k][l] );
          }
          dEdF_tensor[p][q][k][l] = sum;
        }
      }
    }
  }

  //--------------------------------------------------------------------------
  // Precompute dFpmdF_tensor(j, n, k, l)
  // dFpmdF(j,n,k,l) = - sum_{i,m}[ Fpinv(j,i)*dFpdF(i,m,k,l)*Fpinv(m,n) ]
  Real dFpmdF_tensor[3][3][3][3];
  for (unsigned int j = 0; j < 3; j++){
    for (unsigned int n = 0; n < 3; n++){
      for (unsigned int k = 0; k < 3; k++){
        for (unsigned int l = 0; l < 3; l++){
          Real sum = 0.0;
          for (unsigned int i = 0; i < 3; i++){ //summation applies to i, m
            for (unsigned int m = 0; m < 3; m++){
              sum += - Fpinv(j, i) * dFpdF_tensor[i][m][k][l] * Fpinv(m, n);
            }
          }
          dFpmdF_tensor[j][n][k][l] = sum;
        }
      }
    }
  }
  
  //--------------------------------------------------------------------------
  // Compute pk_jacobian using the precomputed tensors.
  RankFourTensor pk_jacobian_val;
  pk_jacobian_val.zero();
  // We sum over the three contributions for every (i,j,k,l):
  // 1. Term1: sum_{m,n} dJpdF(k,l) * _Fe(i,m) * _S(m,n) * Fpinv(j,n)
  // 1. Term1: sum_{m,n} dFedF(i, m, k, l) * _S(m,n) * Fpinv(j,n)
  // 2. Term2: sum_{m,n,p,q} _Fe(i,m)*_C(m,n,p,q)*dEdF(p,q,k,l)*Fpinv(j,n)
  // 3. Term3: sum_{m,n} _Fe(i,m)*_S(m,n)*dFpmdF(j,n,k,l)
  for (unsigned int i = 0; i < 3; i++){
    for (unsigned int j = 0; j < 3; j++){
      for (unsigned int k = 0; k < 3; k++){
        for (unsigned int l = 0; l < 3; l++){
          Real accum = 0.0;
          // Term 1
          for (unsigned int m = 0; m < 3; m++){
            for (unsigned int n = 0; n < 3; n++){
              accum += dJpdF_tensor[k][l] * _Fe[_qp](i, m) * _S[_qp](m, n) * Fpinv(j, n);
            }
          }
          // Term 2
          for (unsigned int m = 0; m < 3; m++){
            for (unsigned int n = 0; n < 3; n++){
              accum += _Jp[_qp] * dFedF_tensor[i][m][k][l] * _S[_qp](m, n) * Fpinv(j, n);
            }
          }
          // Term 3
          for (unsigned int m = 0; m < 3; m++){
            for (unsigned int n = 0; n < 3; n++){
              for (unsigned int p = 0; p < 3; p++){
                for (unsigned int q = 0; q < 3; q++){
                  accum += _Jp[_qp] * _Fe[_qp](i, m) * _C[_qp](m, n, p, q) * dEdF_tensor[p][q][k][l] * Fpinv(j, n);
                }
              }
            }
          }
          // Term 4
          for (unsigned int m = 0; m < 3; m++){
            for (unsigned int n = 0; n < 3; n++){
              accum += _Jp[_qp] * _Fe[_qp](i, m) * _S[_qp](m, n) * dFpmdF_tensor[j][n][k][l];
            }
          }
          pk_jacobian_val(i, j, k, l) = accum;
        }
      }
    }
  }

  //--------------------------------------------------------------------------
  // Off-Diagnoal jacobian: derivatives with respect to pore pressure.
  // --------------------------------------------------------------------------

  // Compute the 2D tensor dFpdP_tensor[i][j] = dFp/dp(i,j)
  // where we approximate: dFp/dp ≈ Fp_dot(i,j) / (p - p_old)
  // Assumes that Fp evolves due to pressure changes (e.g., through Dp evolution)
  // --------------------------------------------------------------------------

  RankTwoTensor dFp_dp; // Tensor to store ∂Fp/∂p
  dFp_dp.zero();        // Initialize to zero
  Real dp = _pore_pressure[_qp] - _pore_pressure_old[_qp]; // Δp
  
  // Approximate derivative: ∂Fp_ij / ∂p ≈ Fp_dot_ij / Δp
  if (_dt != 0.0 && std::abs(dp) > 1e-12){
    for (unsigned int i = 0; i < 3; ++i){
      for (unsigned int j = 0; j < 3; ++j){
        dFp_dp(i,j) = _Fp_dot[_qp](i,j) / dp;
      }
    }
  }
  else{
    dFp_dp.zero();
  }

  // --------------------------------------------------------------------------
  // Compute dJp/dp (derivative of J^p with respect to pore pressure)
  // --------------------------------------------------------------------------

  Real dJp_dp = 0.0;

  // Use chain rule: dJp/dp = Jp * tr(Fpinv^T * dFp/dp)
  if (_dt != 0.0 && std::abs(dp) > 1e-12){
    for (unsigned int i = 0; i < 3; ++i){
      for (unsigned int j = 0; j < 3; ++j){
        dJp_dp += _Jp[_qp] * Fpinv(j, i) * dFp_dp(i, j);
      }
    }
  }
  else{
    dJp_dp = 0.0;
  }

  // --------------------------------------------------------------------------
  // Compute the derivative of effective second Piola-Kirchhoff stress (S) 
  // with respect to pore pressure (pore_pressure)
  // Using: S = S' - α * p
  // So: dS/dp = -α_eff * I
  // --------------------------------------------------------------------------

  RankTwoTensor dS_dp; 
  dS_dp.zero();  // Initialize to zero

  // Compute the effective Biot coefficient (α_eff) based on breakage variable
  const Real B  = _B_breakagevar[_qp];
  const Real alpha_s = _Biot_coeff_s[_qp];
  const Real alpha_g = _Biot_coeff_g[_qp];
  const Real Ms = _Biot_modulus_s[_qp];
  const Real Mg = _Biot_modulus_g[_qp];

  // Effective numerator and denominator for α_eff computation
  Real numerator   = (1.0 - B) * alpha_s * Ms + B * alpha_g * Mg;
  Real denominator = (1.0 - B) * Ms + B * Mg;

  // Final derivative: dS/dp = -α_eff * I
  dS_dp = -(numerator / denominator) * RankTwoTensor::Identity();

  // --------------------------------------------------------------------------
  // Compute dFp_inv/dpore_pressure using the identity:
  // d(Fp⁻¹)/dp = -Fp⁻¹ * dFp/dp * Fp⁻¹
  // --------------------------------------------------------------------------

  RankTwoTensor dFpinv_dp; 
  dFpinv_dp.zero();  // Initialize to zero

  for (unsigned int i = 0; i < 3; ++i){
    for (unsigned int j = 0; j < 3; ++j){
      for (unsigned int k = 0; k < 3; ++k){
        for (unsigned int l = 0; l < 3; ++l){
          dFpinv_dp(i, j) -= Fpinv(i, k) * dFp_dp(k, l) * Fpinv(l, j);
     }
    }
   }
  }

  // --------------------------------------------------------------------------
  // Compute dFe/dpore_pressure using the chain rule:
  // F^e = F · (F^p)⁻¹
  // ⇒ dF^e/dp = F · d(F^p⁻¹)/dp
  // --------------------------------------------------------------------------

  RankTwoTensor dFe_dp;
  dFe_dp.zero();  // initialize

  for (unsigned int i = 0; i < 3; ++i) {
    for (unsigned int j = 0; j < 3; ++j) {
      for (unsigned int k = 0; k < 3; ++k) {
        dFe_dp(i, j) += _F[_qp](i, k) * dFpinv_dp(k, j);
      }
    }
  }

  // --------------------------------------------------------------------------
  // Compute Off-diagonal PK1 Jacobian ∂P/∂pore_pressure at quadrature point
  // --------------------------------------------------------------------------
  RankTwoTensor pk1_off_diag_jacobian_val;
  pk1_off_diag_jacobian_val.zero();

  for (unsigned int i = 0; i < 3; ++i){
    for (unsigned int j = 0; j < 3; ++j){
      Real sum = 0.0;
      for (unsigned int m = 0; m < 3; ++m){
        for (unsigned int n = 0; n < 3; ++n){
          // Contribution 1: dJp/dp * Fe * S * Fp_inv
          Real term1 = dJp_dp * _Fe[_qp](i, m) * _S[_qp](m, n) * Fpinv(n, j);
          // Contribution 2: Jp * dFe/dp * S * Fp_inv
          Real term2 = _Jp[_qp] * dFe_dp(i, m) * _S[_qp](m, n) * Fpinv(n, j);
          // Contribution 3: Jp * Fe * dS/dp * Fp_inv
          Real term3 = _Jp[_qp] * _Fe[_qp](i, m) * dS_dp(m, n) * Fpinv(n, j);
          // Contribution 4: Jp * Fe * S * dFp_inv/dp
          Real term4 = _Jp[_qp] * _Fe[_qp](i, m) * _S[_qp](m, n) * dFpinv_dp(n, j);
          // Sum all contributions
          sum += term1 + term2 + term3 + term4;
        }
      }
      // Store the result
      pk1_off_diag_jacobian_val(i, j) = sum;
    }
  }

  //--------------------------------------------------------------------------
  // PK2-to-PK1 wrapping: using large kinematics formulation.
  if (_large_kinematics)
  {
    // Compute PK1 stress: P = Jp * Fe * S * (Fpinv)^T
    _pk1_stress[_qp] = _Jp[_qp] * _Fe[_qp] * _S[_qp] * Fpinv.transpose();
    _pk1_jacobian[_qp] = pk_jacobian_val;
    _pk1_off_diag_jacobian[_qp] = pk1_off_diag_jacobian_val;
  }
  else
  {
    mooseError("Must select 'large_kinematics' option!");
  }

  // --------------------------------------------------------------------------
  // storing some derivatives for the use of other kernels jacobians
  // --------------------------------------------------------------------------

  // --------------------------------------------------------------------------
  // Compute dI1/dF_{k,l} = sum_{i,j} Fe(i,j) * dFe/dF(i,j,k,l)
  // --------------------------------------------------------------------------

  for (unsigned int k = 0; k < 3; ++k){
    for (unsigned int l = 0; l < 3; ++l){
      Real sum = 0.0;
      for (unsigned int i = 0; i < 3; ++i){
        for (unsigned int j = 0; j < 3; ++j){
          sum += _Fe[_qp](i,j) * dFedF_tensor[i][j][k][l];
        }
      }
      _dI1dF[_qp](k,l) = sum;  
    }
  }

  // --------------------------------------------------------------------------
  // Save dJp/dF tensor to material property _dJpdF[_qp]
  // --------------------------------------------------------------------------
  for (unsigned int k = 0; k < 3; ++k){
    for (unsigned int l = 0; l < 3; ++l){
      _dJpdF[_qp](k, l) = dJpdF_tensor[k][l];
    }
  }

  // --------------------------------------------------------------------------
  // Save dJp/dp scalar to material property _dJpdp[_qp]
  // --------------------------------------------------------------------------
  _dJpdp[_qp] = dJp_dp;

  // --------------------------------------------------------------------------
  // Compute dDpdF_tensor(i,j,k,l) = ∂Dp_ij / ∂F_kl
  // Using: Dp = (I - Fp_inv * Fp_old) / dt
  // ⇒ ∂Dp / ∂F = - (∂Fp_inv/∂F) * Fp_old / dt
  // ⇒ dFp_inv/dF = -Fp⁻¹ * dFp/dF * Fp⁻¹
  // --------------------------------------------------------------------------

  for (unsigned int i = 0; i < 3; ++i){
    for (unsigned int j = 0; j < 3; ++j){
      for (unsigned int k = 0; k < 3; ++k){
        for (unsigned int l = 0; l < 3; ++l){

          Real val = 0.0;
          for (unsigned int m = 0; m < 3; ++m){
            for (unsigned int n = 0; n < 3; ++n){
              for (unsigned int r = 0; r < 3; ++r){
                for (unsigned int s = 0; s < 3; ++s){
                  val += - (1.0 / _dt) * Fpinv(i, m) * dFpdF_tensor[m][n][k][l] * Fpinv(n, r) * _Fp_old[_qp](r, j);
                }
              }
            }
          }

          _dDpdF[_qp](i, j, k, l) = val;
        }
      }
    }
  }

  // --------------------------------------------------------------------------
  // Compute dDpdp_tensor(i,j) = ∂Dp_ij / ∂p
  // Using: Dp = (I - Fp_inv * Fp_old) / dt
  // ⇒ ∂Dp / ∂p = - (∂Fp_inv/∂p) * Fp_old / dt
  // --------------------------------------------------------------------------

  for (unsigned int i = 0; i < 3; ++i){
    for (unsigned int j = 0; j < 3; ++j){
      Real val = 0.0;
      for (unsigned int k = 0; k < 3; ++k){
        val += - dFpinv_dp(i, k) * _Fp_old[_qp](k, j);
      }
      _dDpdp[_qp](i, j) = val / _dt;
    }
  }




}

void
ComputeLagrangianDamageBreakageStressPK2Diffused::computeQpPK2Stress()
{
  /* Evaluate Fp */
  RankTwoTensor Fp_updated = computeQpFp();

  /* Compute Fe */
  RankTwoTensor Fe = _F[_qp] * Fp_updated.inverse();

  /* Compute Ee */
  RankTwoTensor Ee = 0.5 * (Fe.transpose() * Fe - RankTwoTensor::Identity());

  /* Compute Ep */
  RankTwoTensor Ep = 0.5 * (Fp_updated.transpose() * Fp_updated - RankTwoTensor::Identity());

  /* Compute E */
  RankTwoTensor E = Fp_updated.transpose() * Ee * Fp_updated + Ep;

  /* Compute I1 */
  Real I1 = Ee.trace();

  /* Compute I2 */
  Real I2 = 0.0;
  for (unsigned int i = 0; i < 3; ++i){
    for (unsigned int j = 0; j < 3; ++j){
      I2 += Ee(i,j) * Ee(i,j);
    }
  }

  /* Compute xi */
  //here we may need to add small number to avoid singularity
  Real xi = (I1) / (std::sqrt(I2));
  //Catch the nan error in the initial solve
  if (std::isnan(xi)){xi = -std::sqrt(3);}

  Real term11 = (1 - _B_breakagevar[_qp]) * _B_breakagevar[_qp] * _Biot_modulus_s[_qp] * _Biot_modulus_g[_qp] * std::pow((_Biot_coeff_s[_qp] - _Biot_coeff_g[_qp]), 2);
  Real term22 = (1 - _B_breakagevar[_qp]) * _Biot_coeff_s[_qp] * _Biot_modulus_s[_qp] + _B_breakagevar[_qp] * _Biot_coeff_g[_qp] * _Biot_modulus_g[_qp];
  Real term33 = (1 - _B_breakagevar[_qp]) * _Biot_modulus_s[_qp] + _B_breakagevar[_qp] * _Biot_modulus_g[_qp];
  
  /* Compute PK2 stress */
  RankTwoTensor sigma_s = (_lambda_const[_qp] - _damaged_modulus[_qp] / xi) * I1 * RankTwoTensor::Identity() + (2 * _shear_modulus[_qp] - _damaged_modulus[_qp] * xi) * Ee;
  RankTwoTensor sigma_b = (2 * _a2[_qp] + _a1[_qp] / xi + 3 * _a3[_qp] * xi) * I1 * RankTwoTensor::Identity() + (2 * _a0[_qp] + _a1[_qp] * xi - _a3[_qp] * std::pow(xi, 3)) * Ee;
  RankTwoTensor fluid_contribution = term11 / term33 * I1 * RankTwoTensor::Identity() - term22 / term33 * _pore_pressure[_qp] * RankTwoTensor::Identity();
  RankTwoTensor sigma_total;
  
  sigma_total = (1 - _B_breakagevar[_qp]) * sigma_s + _B_breakagevar[_qp] * sigma_b + fluid_contribution;

  //save
  _Ep[_qp] = Ep;
  _S[_qp] = sigma_total;

  /* Compute plastic stress */
  _Tp[_qp] = Fe.transpose() * Fe * sigma_total;

  /* Compute the effective plastic stress from thermodynamics */
  RankTwoTensor Tp_eff;
  Real Jp = _Fp[_qp].det();
  Tp_eff = _Tp[_qp] + Jp * _pore_pressure[_qp] * RankTwoTensor::Identity();

  /* Compute deviatoric stress tensor */
  _Tau[_qp] = Tp_eff - 0.3333 * (Tp_eff.trace()) * RankTwoTensor::Identity();

  /* Compute tangent */
  RankFourTensor tangent;
  computeQpTangentModulus(tangent,I1,I2,xi,Ee);

  //save
  _C[_qp] = tangent;

  /* Save other parameters */
  _Fp[_qp] = Fp_updated;
  _Fe[_qp] = Fe;
  _Ee[_qp] = Ee;
  _E[_qp]  = E;
  _I1[_qp] = I1;
  _I2[_qp] = I2;
  _xi[_qp] = xi;

  //compute deviatoric strain rate
  computeDeviatroicStrainRateTensor();

}

RankTwoTensor
ComputeLagrangianDamageBreakageStressPK2Diffused::computeQpFp()
{
  // //NOT FINISHED
  // //Compute eigen-decomposition of Tau
  // std::vector<Real> eigval(3, 0.0);
  // RankTwoTensor diag;
  // RankTwoTensor Q;
  // RankTwoTensor PowTau;

  // _Tau_old[_qp].symmetricEigenvaluesEigenvectors(eigval, Q);

  // const Real eps = 1e-8;
  // for (unsigned int i = 0; i < 3; ++i)
  // {
  //   Real eig = std::max(std::abs(eigval[i]), eps); // clamp
  //   diag(i, i) = std::copysign(std::pow(eig, _m2[_qp]), eigval[i]); // preserve sign
  // }

  // PowTau = Q * diag * Q.transpose();

  // //Get old Tau
  // RankTwoTensor Tau_old = PowTau;

  // //Get equvialent deviatroic stress scalar
  // Real Tau_eq = 0.0;
  // for (unsigned int p = 0; p < 3; p++){
  //   for (unsigned int q = 0; q < 3; q++){
  //     Tau_eq += 2.0/3.0 * Tau_old(p,q) * Tau_old(p,q);
  //   }
  // }

  // Tau_eq = std::sqrt(Tau_eq); 

  // //Get deviatroic stress direction
  // RankTwoTensor N; N.zero();
  // // Epsilon to avoid division by zero
  // if (Tau_eq != 0.0){
  //   //Compute deviatroic stress direction
  //   for (unsigned int p = 0; p < 3; p++){
  //     for (unsigned int q = 0; q < 3; q++){
  //       N(p,q) = Tau_old(p,q) / Tau_eq;
  //     }
  //   }
  // }

  //Apply power operation on every element of Tau
  //let's assume m2 = 1, and not apply pow on its elements
  RankTwoTensor Tau_old_power_m2 = _Tau_old[_qp];

  // Define equivalent shear rate nu
  if (_use_state_var_evolution_mat[_qp])
  {
    computestatedependentDp();
  }
  else if (_use_dilatancy)
  {
    _shear_rate_nu[_qp] = _C_g[_qp] * std::pow(_B_breakagevar_old[_qp], _m1[_qp]) * Tau_old_power_m2;

    _eta[_qp] = _eta_old[_qp] + _dilatancy_function_beta[_qp] * _C_g[_qp] * std::pow(_B_breakagevar_old[_qp], _m1[_qp]) * _dt;
    
    _dilatancy_function_beta[_qp] = _anand_param_go_mat * std::pow( 1 - _eta[_qp] / _anand_param_eta_cv_mat, _anand_param_p_mat );   
  }
  else
  {
    _shear_rate_nu[_qp] = _C_g[_qp] * std::pow(_B_breakagevar_old[_qp], _m1[_qp]) * Tau_old_power_m2;
  }

  //Compute Plastic Deformation Rate Tensor Dp at t_{n+1} using quantities from t_{n}
  RankTwoTensor Dp = _shear_rate_nu[_qp]; 

  if (_use_dilatancy)
  {
    Dp += _C_g[_qp] * std::pow(_B_breakagevar_old[_qp], _m1[_qp]) * _dilatancy_function_beta[_qp] / 3 * RankTwoTensor::Identity(); 
  }
 
  //Compute Cp = I - Dp dt
  RankTwoTensor Cp = RankTwoTensor::Identity() - Dp * _dt;

  //Use Implicit Euler Integration, Update Fp
  RankTwoTensor Fp_updated = Cp.inverse() * _Fp_old[_qp];

  //Save Plastic Deformation Rate Tensor
  _Dp[_qp] = Dp;

  return Fp_updated;
}

void
ComputeLagrangianDamageBreakageStressPK2Diffused::computeQpTangentModulus(RankFourTensor & tangent, 
                                                                  Real I1, 
                                                                  Real I2, 
                                                                  Real xi, 
                                                                  RankTwoTensor Ee)
{

  // Use consistent values - same as in stress computation

  // Use the SAME values as in stress computation
  Real lambda_out = _lambda_const[_qp];
  Real shear_modulus_out = _shear_modulus[_qp];
  Real gamma_damaged_out = _damaged_modulus[_qp];

  Real a0 = _a0[_qp];
  Real a1 = _a1[_qp];
  Real a2 = _a2[_qp];
  Real a3 = _a3[_qp];

  Real term11 = (1 - _B_breakagevar[_qp]) * _B_breakagevar[_qp] * _Biot_modulus_s[_qp] * _Biot_modulus_g[_qp] * std::pow((_Biot_coeff_s[_qp] - _Biot_coeff_g[_qp]), 2);
  Real term33 = (1 - _B_breakagevar[_qp]) * _Biot_modulus_s[_qp] + _B_breakagevar[_qp] * _Biot_modulus_g[_qp];
  
  // Safety check for small I2
  const Real adjusted_I2 = std::max(I2, 1e-12);
  const Real sqrt_I2 = std::sqrt(adjusted_I2);
  const RankTwoTensor identity = RankTwoTensor::Identity();

  // Check for limiting case: alpha = 0, B = 0 (should return elasticity tensor)
  // if (std::abs(_alpha_damagedvar_aux[_qp]) < 1e-12 && std::abs(_B_damagedvar_aux[_qp]) < 1e-12) {
  //   // Standard elasticity tensor: C_ijkl = λ δ_ij δ_kl + μ (δ_ik δ_jl + δ_il δ_jk)
  //   tangent.zero();
  //   for (unsigned int i = 0; i < 3; ++i) {
  //     for (unsigned int j = 0; j < 3; ++j) {
  //       for (unsigned int k = 0; k < 3; ++k) {
  //         for (unsigned int l = 0; l < 3; ++l) {
  //           tangent(i, j, k, l) = _lambda_o * identity(i, j) * identity(k, l) + 
  //                                 _shear_modulus_o * (identity(i, k) * identity(j, l) + identity(i, l) * identity(j, k));
  //         }
  //       }
  //     }
  //   }
  //   return;
  // }

  // Corrected derivative: ∂ξ/∂E_kl = ∂(I1/√I2)/∂E_kl
  // = (∂I1/∂E_kl * √I2 - I1 * ∂I2/∂E_kl / (2√I2)) / I2
  // where ∂I1/∂E_kl = δ_kl and ∂I2/∂E_kl = 2*E_kl
  RankTwoTensor dxidE_tensor;
  for (unsigned int k = 0; k < 3; ++k) {
    for (unsigned int l = 0; l < 3; ++l) {
      dxidE_tensor(k, l) = identity(k, l) / sqrt_I2 - I1 * Ee(k, l) / std::pow(adjusted_I2, 1.5);
    }
  }

  // ∂(1/ξ)/∂E = -1/ξ² * ∂ξ/∂E
  const RankTwoTensor dxim1dE_tensor = dxidE_tensor * (-1.0 / (xi * xi));

  // Compute solid phase tangent (dSs/dE)
  const Real lambda_term = lambda_out - gamma_damaged_out / xi;
  const Real shear_term = 2.0 * shear_modulus_out - gamma_damaged_out * xi;

  RankFourTensor dSsdE;
  dSsdE.zero();
  
  // CORRECTED: Complete implementation of solid phase tangent
  // ∂S^s_ij/∂E_kl = (-γ ∂ξ^(-1)/∂E_kl)I_1 δ_ij + (λ - γ/ξ) ∂I_1/∂E_kl δ_ij + (-γ ∂ξ/∂E_kl)E_ij + (2μ - γξ) ∂E_ij/∂E_kl
  for (unsigned int i = 0; i < 3; ++i) {
    for (unsigned int j = 0; j < 3; ++j) {
      for (unsigned int k = 0; k < 3; ++k) {
        for (unsigned int l = 0; l < 3; ++l) {
          // Term 1a: (λ - γ/ξ) * ∂I1/∂E_kl * δ_ij = (λ - γ/ξ) * δ_kl * δ_ij
          dSsdE(i, j, k, l) += lambda_term * identity(i, j) * identity(k, l);
          
          // Term 1b: (-γ ∂ξ^(-1)/∂E_kl) * I1 * δ_ij - PREVIOUSLY MISSING
          dSsdE(i, j, k, l) -= gamma_damaged_out * dxim1dE_tensor(k, l) * I1 * identity(i, j);
          
          // Term 2a: (2μ - γξ) * ∂E_ij/∂E_kl
          Real I4_ijkl = 0.5 * (identity(i, k) * identity(j, l) + identity(i, l) * identity(j, k));
          dSsdE(i, j, k, l) += shear_term * I4_ijkl;
          
          // Term 2b: (-γ ∂ξ/∂E_kl) * E_ij
          dSsdE(i, j, k, l) -= gamma_damaged_out * dxidE_tensor(k, l) * Ee(i, j);
        }
      }
    }
  }

  // Compute granular phase tangent (dSb/dE)
  const Real coeff2_b = 2.0 * a2 + a1 / xi + 3.0 * a3 * xi;
  const Real coeff4_b = 2.0 * a0 + a1 * xi - a3 * xi * xi * xi;

  RankFourTensor dSbdE;
  dSbdE.zero();
  
 // CORRECTED: Complete implementation of granular phase tangent
  // ∂S^b_ij/∂E_kl = (a_1 ∂ξ^(-1)/∂E_kl + 3a_3 ∂ξ/∂E_kl)I_1 δ_ij + (2a_2 + a_1/ξ + 3a_3ξ) ∂I_1/∂E_kl δ_ij
  //                + (a_1 ∂ξ/∂E_kl - a_3 ∂ξ^3/∂E_kl)E_ij + (2a_0 + a_1ξ - a_3ξ^3) ∂E_ij/∂E_kl
  for (unsigned int i = 0; i < 3; ++i) {
    for (unsigned int j = 0; j < 3; ++j) {
      for (unsigned int k = 0; k < 3; ++k) {
        for (unsigned int l = 0; l < 3; ++l) {
          // Term 1a: (2a_2 + a_1/ξ + 3a_3ξ) * ∂I1/∂E_kl * δ_ij = coeff2_b * δ_kl * δ_ij
          dSbdE(i, j, k, l) += coeff2_b * identity(i, j) * identity(k, l);
          
          // Term 1b: a_1 * ∂ξ^(-1)/∂E_kl * I1 * δ_ij - PREVIOUSLY MISSING
          dSbdE(i, j, k, l) += a1 * dxim1dE_tensor(k, l) * I1 * identity(i, j);
          
          // Term 1c: 3a_3 * ∂ξ/∂E_kl * I1 * δ_ij
          dSbdE(i, j, k, l) += 3.0 * a3 * dxidE_tensor(k, l) * I1 * identity(i, j);
          
          // Term 2a: (2a_0 + a_1ξ - a_3ξ^3) * ∂E_ij/∂E_kl
          Real I4_ijkl = 0.5 * (identity(i, k) * identity(j, l) + identity(i, l) * identity(j, k));
          dSbdE(i, j, k, l) += coeff4_b * I4_ijkl;
          
          // Term 2b: a_1 * ∂ξ/∂E_kl * E_ij
          dSbdE(i, j, k, l) += a1 * dxidE_tensor(k, l) * Ee(i, j);
          
          // Term 2c: -a_3 * ∂ξ^3/∂E_kl * E_ij
          // ∂ξ^3/∂E_kl = 3ξ^2 * ∂ξ/∂E_kl
          dSbdE(i, j, k, l) -= a3 * 3.0 * xi * xi * dxidE_tensor(k, l) * Ee(i, j);
        }
      }
    }
  }

  // Combine: tangent = (1-B)*dSs/dE + B*dSb/dE
  tangent = dSsdE * (1.0 - _B_breakagevar[_qp]) + dSbdE * _B_breakagevar[_qp]; 

  // Compute fluid tangent contribution: d(fluid_contribution)/dE
  // fluid_contribution = term11/term33 * I1 * Identity - term22/term33 * p * Identity
  // d(fluid_contribution)/dE = term11/term33 * dI1/dE * Identity
  // where dI1/dE = Identity (trace operation)
    
  Real fluid_contrib = term11 / term33;
    
  RankFourTensor fluid_tangent;
  fluid_tangent.zero();
    
  // d(fluid_contribution_ij)/dE_kl = fluid_coeff * d(I1)/dE_kl * δ_ij
  // where d(I1)/dE_kl = δ_kl
  for (unsigned int i = 0; i < 3; ++i) {
    for (unsigned int j = 0; j < 3; ++j) {
      for (unsigned int k = 0; k < 3; ++k) {
        for (unsigned int l = 0; l < 3; ++l) {
            fluid_tangent(i, j, k, l) = fluid_contrib * identity(k, l) * identity(i, j);
        }
      }
    }
  }
    
  // Add fluid tangent to total tangent
  tangent += fluid_tangent;

}

void 
ComputeLagrangianDamageBreakageStressPK2Diffused::computeDmatrix()
{
  //Compute deformation rate D
  _D[_qp] = 0.5 * ( _velgrad_L[_qp] + _velgrad_L[_qp].transpose() );
}

void
ComputeLagrangianDamageBreakageStressPK2Diffused::computeDeviatroicStrainRateTensor()
{
  //Compute D matrix
  computeDmatrix();
  //Compute strain rate E_dot = F^T * D * F
  RankTwoTensor E_dot = _F[_qp].transpose() * _D[_qp] * _F[_qp];
  //Compute deviatoric strain rate tensor E_dev_dot 
  RankTwoTensor E_dev_dot = E_dot - (1.0/3.0) * E_dot.trace() * RankTwoTensor::Identity();
  //Compute J2_dot = 1/2 * E_dev_dot(i,j) * E_dev_dot(i,j)
  Real J2_dot = 0.0;
  for (unsigned int i = 0; i < 3; ++i){
    for (unsigned int j = 0; j < 3; ++j){
      J2_dot += 0.5 * E_dev_dot(i,j) * E_dev_dot(i,j);
    }
  }
  //Compute equivalent strain rate
  _deviatroic_strain_rate[_qp] = std::sqrt(2.0/3.0 * J2_dot);
}

void 
ComputeLagrangianDamageBreakageStressPK2Diffused::computestatedependentDp()
{
  //Compute equivalent plastic strain rate
  //------------------------------------------------------------------------------//
  Real Dp_eq = 0.0;
  for (unsigned int i = 0; i < 3; i++){
    for (unsigned int j = 0; j < 3; j++){
      Dp_eq += _Dp[_qp](i,j) * _Dp[_qp](i,j);
    }
  }

  Dp_eq = std::sqrt(2.0/3.0 * Dp_eq);
  //------------------------------------------------------------------------------//

  //NOT FINISHED
  //------------------------------------------------------------------------------//
  //Compute eigen-decomposition of Tau
  std::vector<Real> eigval(3, 0.0);
  RankTwoTensor diag;
  RankTwoTensor Q;
  RankTwoTensor PowTau;

  _Tau_old[_qp].symmetricEigenvaluesEigenvectors(eigval, Q);

  const Real eps = 1e-8;
  for (unsigned int i = 0; i < 3; ++i)
  {
    Real eig = std::max(std::abs(eigval[i]), eps); // clamp
    diag(i, i) = std::copysign(std::pow(eig, _m2[_qp]), eigval[i]); // preserve sign
  }

  PowTau = Q * diag * Q.transpose();

  //Get old Tau
  RankTwoTensor Tau_old = PowTau;

  //Get equvialent deviatroic stress scalar
  Real Tau_eq = 0.0;
  for (unsigned int p = 0; p < 3; p++){
    for (unsigned int q = 0; q < 3; q++){
      Tau_eq += 2.0/3.0 * Tau_old(p,q) * Tau_old(p,q);
    }
  }

  Tau_eq = std::sqrt(Tau_eq); 
  //------------------------------------------------------------------------------//

  //Compute shear rate nu
  _Theta[_qp] = _dt * ( 1 - 4 * Dp_eq * _Theta_old[_qp] ) + _Theta_old[_qp];
  _shear_rate_nu[_qp] = _C_g[_qp] * std::pow(_B_breakagevar_old[_qp],_m1[_qp]) * std::pow(Tau_eq,_m2[_qp]) * std::pow((1.0*_Theta[_qp])/(1.0*_const_theta_o_mat[_qp]),(-1.0*_const_B_mat[_qp])/(1.0*_const_A_mat[_qp]));
}

