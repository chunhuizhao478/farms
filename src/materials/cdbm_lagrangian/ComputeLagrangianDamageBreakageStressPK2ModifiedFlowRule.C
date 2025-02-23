//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ComputeLagrangianDamageBreakageStressPK2ModifiedFlowRule.h"

/// Native interface for providing the 2nd Piola Kirchhoff stress
///
/// This class *implements* the 2nd PK stress update, providing:
///   1) The 2nd PK stress
///   2) The derivative of the 2nd PK stress wrt the Cauchy-Green strain
///
/// and wraps these to provide:
///   1) The 1st PK stress
///   2) d(PK1)/d(F)
///
/// Modified Flow Rule, Created By Chunhui Zhao, Feb 22th 2025
///
/// Introduce state variable theta to control a reasonable residual stress level
/// 
/// (1) Define material property state variable \theta
///
/// - Handle new constant material property passed from DamageBreakageMaterial: const_A const_B const_theta_o 
///
/// (2) The plastic strain rate:
/// D^p_{ij} = Cg B^{m_1} exp(T^p_{o,ij}/A) (\theta/\theta_o)^{-B/A}
///
/// - Compute exponential of real symmetric tensor T^p_o using spectral decomposition: exp(T^p_{o,ij}/A) = Q exp(diag(T^p_{o,ij}/A)) Q^T
///
/// Or using without exponential form: D^p_{ij} = Cg B^{m_1} T^p_{o,ij} (\theta/\theta_o)^{-B/A}
///
/// (3) The evolution of state variable:
/// dot{theta} = 1 - D^p_{eq} theta
///
/// - Compute using implicit time integration: theta^{t+1} = (theta^t + \Delta t)/(1 + D^p_eq \Delta t)
/// - Compute scalar quantity D^p_{eq} = sqrt(2/3*D^p_{ij}*D^p_{ij})
///

registerMooseObject("farmsApp", ComputeLagrangianDamageBreakageStressPK2ModifiedFlowRule);

InputParameters
ComputeLagrangianDamageBreakageStressPK2ModifiedFlowRule::validParams()
{
  InputParameters params = ComputeLagrangianStressPK1::validParams();
  return params;
}

ComputeLagrangianDamageBreakageStressPK2ModifiedFlowRule::ComputeLagrangianDamageBreakageStressPK2ModifiedFlowRule(const InputParameters & parameters)
  : ComputeLagrangianStressPK1(parameters),
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
  _Theta(declareProperty<Real>(_base_name + "state_variable")),
  _lambda_const(getMaterialProperty<Real>("lambda_const")),
  _shear_modulus(getMaterialProperty<Real>("shear_modulus")),
  _damaged_modulus(getMaterialProperty<Real>("damaged_modulus")),
  _B_breakagevar(getMaterialProperty<Real>("B_damagedvar")),
  _B_breakagevar_old(getMaterialPropertyOldByName<Real>("B_damagedvar")),
  _Tau_old(getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "deviatroic_stress")),
  _Fp_old(getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "plastic_deformation_gradient")),
  _F_old(getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "deformation_gradient")),
  _Ep_old(getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "plastic_strain")),
  _Theta_old(getMaterialPropertyOldByName<Real>(_base_name + "state_variable")),
  _Dp_old(getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "plastic_strain_rate")),
  _C_g(getMaterialProperty<Real>("C_g")),
  _m1(getMaterialProperty<Real>("m1")),
  _m2(getMaterialProperty<Real>("m2")),
  _a0(getMaterialProperty<Real>("a0")),
  _a1(getMaterialProperty<Real>("a1")),
  _a2(getMaterialProperty<Real>("a2")),
  _a3(getMaterialProperty<Real>("a3")),
  _use_state_var_evolution_mat(getMaterialProperty<bool>("use_state_var_evolution_mat")),
  _const_A_mat(getMaterialProperty<Real>("const_A_mat")),
  _const_B_mat(getMaterialProperty<Real>("const_B_mat")),
  _const_theta_o_mat(getMaterialProperty<Real>("const_theta_o_mat")),
  _use_vels_build_L_mat(getMaterialProperty<bool>("use_vels_build_L_mat")),
  _velgrad_L(getMaterialProperty<RankTwoTensor>("velgrad_L"))
  //_dim(_mesh.dimension())
{
}

//Rules:See https://github.com/idaholab/moose/discussions/19450
//Only the object that declares the material property can assign values to it.
//Objects can request material properties, gaining read-only access to their values.
//When any object (including the object that declares it) requests the old value of a material property, that property becomes "stateful".
//All stateful material properties must be initialized within the initQpStatefulProperties call. 
//
void
ComputeLagrangianDamageBreakageStressPK2ModifiedFlowRule::initQpStatefulProperties()
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
  _Theta[_qp] = _const_theta_o_mat[_qp];
  _Dp[_qp].zero();
}

void
ComputeLagrangianDamageBreakageStressPK2ModifiedFlowRule::computeQpPK1Stress()
{
  // PK2 update
  computeQpPK2Stress();

  // Compute Jp, Fp^{-1} //depends on the plastic deformation rate, here no volumetric strain upodate, Jp = 1 (checked)
  _Jp[_qp] = _Fp[_qp].det();
  RankTwoTensor Fpinv = _Fp[_qp].inverse();

  //Compute deformation rate D
  _D[_qp] = 0.5 * ( _velgrad_L[_qp] + _velgrad_L[_qp].transpose() );

  // Compute Fp_dot, F_dot
  // Here we approximate the rate by first-order, not sure if this is sufficient for varying time steps
  // currently MOOSE don't support getMaterialPropertyDot
  RankTwoTensor Fp_dot;
  RankTwoTensor F_dot;
  
  if (_use_vels_build_L_mat[_qp]){ //use true deformation rate

    for (unsigned int i = 0; i < 3; i++){
      for (unsigned int j = 0; j < 3; j++){
        for (unsigned int m = 0; m < 3; m++){
          Fp_dot(i,j) += _Dp[_qp](i,m) * _Fp[_qp](m,j);
          F_dot(i,j) += _D[_qp](i,m) * _F[_qp](m,j); 
        }
      }
    } 

  }
  else{ //finite difference approximation

    for (unsigned int i = 0; i < 3; i++){
      for (unsigned int j = 0; j < 3; j++){
          F_dot(i,j)  = (_F[_qp](i,j) - _F_old[_qp](i,j) ) / _dt; 
          //F_dot(i,j)  = (_F[_qp](i,j) - _F_old[_qp](i,j) ); 
          for (unsigned int m = 0; m < 3; m++){
            Fp_dot(i,j) += _Dp[_qp](i,m) * _Fp[_qp](m,j);
          }
      }
    }
  
  }
 
  // Compute delta function
  auto delta = [](int i, int j) -> Real {
    return (i == j) ? 1.0 : 0.0;
  };  

  //Compute dFpdF //need to confirm
  auto dFpdF = [&](int i, int j, int k, int l) -> Real {
    if (_dt == 0.0){ //Steady, set zero
      return 0.0;
    } 
    else{ //Transient
      if (Fp_dot(i,j) == 0.0 || F_dot(k,l) == 0.0 ){ //no change of viscoelastic dg, set zero
        return 0.0;
      }
      else{
        return Fp_dot(i,j)/F_dot(k,l);
      }
    }
  };

  // Compute dFedF
  auto dFedF = [&](int i, int m, int k, int l) -> Real {

    // Since dFp/dF is zero (Fp is explicit), dFedF simplifies
    Real dFedF_val = delta(i,k) * Fpinv(l,m);

    //here apply summations to {h,r}
    for (unsigned int h = 0; h < 3; h++){
      for (unsigned int r = 0; r < 3; r++){
        dFedF_val -= _Fe[_qp](i,h) * dFpdF(h,r,k,l) * Fpinv(r,m);
      }
    }    

    return dFedF_val;

  }; 

  // Compute dEdF
  auto dEdF = [&](int p, int q, int k, int l) -> Real {
    
    //initialize value
    Real dEdF_val = 0.0;

    //here apply summations to {m}
    for (unsigned int m = 0; m < 3; m++){
      dEdF_val += 0.5 * ( dFedF(m,p,k,l) * _Fe[_qp](m,q) + _Fe[_qp](m,p) * dFedF(m,q,k,l) );
    }

    return dEdF_val;

  };

  // Compute dFpmdF
  auto dFpmdF = [&](int j, int n, int k, int l) -> Real {

    //initialize value
    Real dFpmdF_val = 0.0;

    //here apply summations to {i,a}
    for (unsigned int i = 0; i < 3; i++){
      for (unsigned int m = 0; m < 3; m++){
        dFpmdF_val += -1.0 * Fpinv(j,i) * dFpdF(i,m,k,l) * Fpinv(m,n);
      }
    }

    return dFpmdF_val;

  };

  //Compute pk_jacobian
  RankFourTensor pk_jacobian_val;
  pk_jacobian_val.zero();  // Make sure the tensor starts with zero values
  for (unsigned int i = 0; i < 3; i++){
    for (unsigned int j = 0; j < 3; j++){
      for (unsigned int k = 0; k < 3; k++){
        for (unsigned int l = 0; l < 3; l++){
          for (unsigned int m = 0; m < 3; m++)
          {
            // First term: dFedF(i,m,k,l) * S(m,n) * Fpinv(j,n)
            for (unsigned int n = 0; n < 3; n++){
              pk_jacobian_val(i,j,k,l) += dFedF(i,m,k,l) * _S[_qp](m,n) * Fpinv(j,n);
            }
            
            // Second term: Fe(i,m) * C(m,n,p,q) * dEdF(p,q,k,l) * Fpinv(j,n)
            for (unsigned int n = 0; n < 3; n++){
              for (unsigned int p = 0; p < 3; p++){
                for (unsigned int q = 0; q < 3; q++){
                  pk_jacobian_val(i,j,k,l) += _Fe[_qp](i,m) * _C[_qp](m,n,p,q) * dEdF(p,q,k,l) * Fpinv(j,n);
                }
              }
            }
            
            // Third term: Fe(i,m) * S(m,n) * dFpmdF(j,n,k,l)
            for (unsigned int n = 0; n < 3; n++){
              pk_jacobian_val(i,j,k,l) += _Fe[_qp](i,m) * _S[_qp](m,n) * dFpmdF(j,n,k,l);
            }
          }
        }
      }
    }
  }

  // Complicated wrapping from PK2 to PK1, see documentation on overleaf
  if (_large_kinematics)
  {
    //if there is plastic
    // Compute pk1 stress
    _pk1_stress[_qp] = _Fe[_qp] * _S[_qp] * Fpinv.transpose();

    // Compute pk1 jacobian
    _pk1_jacobian[_qp] = pk_jacobian_val;

  }
  else
  {
    mooseError("Must selection 'large_kinematics' option!");
  }

}


void
ComputeLagrangianDamageBreakageStressPK2ModifiedFlowRule::computeQpPK2Stress()
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

  /* Compute PK2 stress */
  RankTwoTensor sigma_s = (_lambda_const[_qp] - _damaged_modulus[_qp] / xi) * I1 * RankTwoTensor::Identity() + (2 * _shear_modulus[_qp] - _damaged_modulus[_qp] * xi) * Ee;
  RankTwoTensor sigma_b = (2 * _a2[_qp] + _a1[_qp] / xi + 3 * _a3[_qp] * xi) * I1 * RankTwoTensor::Identity() + (2 * _a0[_qp] + _a1[_qp] * xi - _a3[_qp] * std::pow(xi, 3)) * Ee;
  RankTwoTensor sigma_total = (1 - _B_breakagevar[_qp]) * sigma_s + _B_breakagevar[_qp] * sigma_b;

  //save
  _Ep[_qp] = Ep;
  _S[_qp] = sigma_total;

  /* Compute plastic stress */
  _Tp[_qp] = Fe.transpose() * Fe * sigma_total;

  //compute deviatroic stress tensor //save
  _Tau[_qp] = _Tp[_qp] - 0.3333 * ( _Tp[_qp].trace() ) * RankTwoTensor::Identity();

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

}

void
ComputeLagrangianDamageBreakageStressPK2ModifiedFlowRule::computeTheta()
{

  //Compute equivalent plastic strain rate
  Real Dp_eq = 0.0;
  for (unsigned int i = 0; i < 3; i++){
    for (unsigned int j = 0; j < 3; j++){
      Dp_eq += _Dp_old[_qp](i,j) * _Dp_old[_qp](i,j);
    }
  }

  Dp_eq = std::sqrt(2.0/3.0 * Dp_eq);

  //Compute theta using explicit time integration
  _Theta[_qp] = _Theta_old[_qp] + (1.0 - Dp_eq * _Theta_old[_qp]) * _dt;
}

RankTwoTensor
ComputeLagrangianDamageBreakageStressPK2ModifiedFlowRule::computeQpFp()
{
  //Compute eigen-decomposition of Tau
  // std::vector<Real> eigval(3, 0.0);
  // RankTwoTensor<Real> diag;
  // RankTwoTensor<Real> Q;
  // RankTwoTensor<Real> ExpTau;

  // _Tau_old[_qp].symmetricEigenvaluesEigenvectors(eigval, Q);

  // diag[0][0] = std::exp(eigval[0]/_const_A_mat[_qp]);
  // diag[1][1] = std::exp(eigval[1]/_const_A_mat[_qp]);
  // diag[2][2] = std::exp(eigval[2]/_const_A_mat[_qp]);

  // ExpTau = Q * diag * Q.transpose();

  //Compute state variable theta
  computeTheta();

  //Compute Plastic Deformation Rate Tensor Dp at t_{n+1} using quantities from t_{n}
  RankTwoTensor Dp;
  if (_use_state_var_evolution_mat[_qp]){
    _Dp[_qp] = _C_g[_qp] * std::pow(_B_breakagevar_old[_qp],_m1[_qp]) * _Tau_old[_qp] * std::pow(_Theta[_qp]/_const_theta_o_mat[_qp],-1.0*_const_B_mat[_qp]/_const_A_mat[_qp]);
  }
  else{
    mooseError("Must select 'use_state_var_evolution_mat = true' option in DamageBreakageMaterial to use this material object");
  }

  //Compute Cp = I - Dp dt
  RankTwoTensor Cp = RankTwoTensor::Identity() - _Dp[_qp] * _dt;

  //Use Implicit Euler Integration, Update Fp
  RankTwoTensor Fp_updated = Cp.inverse() * _Fp_old[_qp];

  return Fp_updated;
}

void
ComputeLagrangianDamageBreakageStressPK2ModifiedFlowRule::computeQpTangentModulus(RankFourTensor & tangent, 
                                                                  Real I1, 
                                                                  Real I2, 
                                                                  Real xi, 
                                                                  RankTwoTensor Ee)
{

  //define functions for derivatives

  //delta function
  auto delta = [](int i, int j) -> Real {
    return (i == j) ? 1.0 : 0.0;
  };

  //dI1_dE_{kl}
  auto dI1dE = [&](int k, int l) -> Real {
    return delta(k,l);
  };

  //dI2_dE_{kl}
  auto dI2dE = [&](int k, int l) -> Real {
    return 2 * Ee(k,l);
  };

  //dxi_dE_{kl}
  auto dxidE = [&](int k, int l) -> Real {
    
    // Epsilon to avoid division by zero
    const Real epsilon = 1e-12;
    // Adjust I2 if necessary
    Real adjusted_I2 = I2;
    if (I2 <= epsilon) {
      //mooseWarning("I2 is zero or too small (I2 = ", I2, "), adjusting to epsilon.");
      adjusted_I2 = epsilon;
    }

    Real dxidE = 0.5 * pow(adjusted_I2,-1.5) * dI2dE(k,l) * I1;
    //mooseInfo("I1 = ", I1, ", I2 = ", I2);
    if (std::isnan(dxidE)){mooseError("dxidE");}
    return delta(k,l) * pow(adjusted_I2,-0.5) - 0.5 * pow(adjusted_I2,-1.5) * dI2dE(k,l) * I1;
  };

  //dE_{ij}_dE_{kl}
  auto dEdE = [&](int i, int j, int k, int l) -> Real {
    return delta(i,k) * delta(j,l);
    //return 0.5 * ( delta(i,k) * delta(j,l) + delta(i,l) * delta(j,k) ); //its symmetric form
  };

  //dxi^{-1}_dE_{kl}
  auto dxim1dE = [&](int k, int l) -> Real {
    return -1.0 * pow(xi,-2.0) * dxidE(k,l);
  };

  //dxi^3_dE_{kl}
  auto dxi3dE = [&](int k, int l) -> Real {
    return 3 * pow(xi,2) * dxidE(k,l);
  };

  //dSe_{ij}_dE_{kl}
  auto dSedE = [&](int i, int j, int k, int l) -> Real {
    Real dSedE_components = (- _damaged_modulus[_qp] * dxim1dE(k,l) ) * I1 * delta(i,j);
    dSedE_components += ( _lambda_const[_qp] - _damaged_modulus[_qp] / xi ) * dI1dE(k,l) * delta(i,j);
    dSedE_components += (- _damaged_modulus[_qp] * dxidE(k,l) ) * Ee(i,j);
    dSedE_components += ( 2 * _shear_modulus[_qp] - _damaged_modulus[_qp] * xi ) * dEdE(i,j,k,l);
    if (std::isnan(dSedE_components)){mooseError("dSedE_components");}
    return dSedE_components;
  };

  //dSb_{ij}_dE_{kl}
  auto dSbdE = [&](int i, int j, int k, int l) -> Real {
    Real dSbdE_components = ( _a1[_qp] * dxim1dE(k,l) + 3 * _a3[_qp] * dxidE(k,l) ) * I1 * delta(i,j);
    dSbdE_components += ( 2 * _a2[_qp] + _a1[_qp] / xi + 3 * _a3[_qp] * xi ) * dI1dE(k,l) * delta(i,j);
    dSbdE_components += ( _a1[_qp] * dxidE(k,l) - _a3[_qp] * dxi3dE(k,l) ) * Ee(i,j);
    dSbdE_components += ( 2 * _a0[_qp] + _a1[_qp] * xi - _a3[_qp] * pow(xi,3) ) * dEdE(i,j,k,l);
    return dSbdE_components;
  };

  //dS_{ij}_dE_{kl}
  auto dSdE = [&](int i, int j, int k, int l) -> Real {
    return (1 - _B_breakagevar[_qp]) * dSedE(i,j,k,l) + _B_breakagevar[_qp] * dSbdE(i,j,k,l);
  };

  // Compute tangent modulus C
  for (unsigned int i = 0; i < 3; i++){
    for (unsigned int j = 0; j < 3; j++){
      for (unsigned int k = 0; k < 3; k++){
        for (unsigned int l = 0; l < 3; l++){
          if (std::isnan(dSdE(i,j,k,l))){mooseError("encounter nan error: dSdE(i,j,k,l)");}
          tangent(i,j,k,l) += dSdE(i,j,k,l);
        }
      }
    }
  }

  // const Real adjusted_I2 = (I2 <= 1e-12) ? 1e-12 : I2;
  // const RankTwoTensor identity = RankTwoTensor::Identity();

  // // Precompute dxidE tensor
  // RankTwoTensor dxidE_tensor;
  // for (unsigned int k = 0; k < 3; ++k)
  //   for (unsigned int l = 0; l < 3; ++l)
  //     dxidE_tensor(k, l) = (identity(k, l) * adjusted_I2 - I1 * Ee(k, l)) / std::pow(adjusted_I2, 1.5);

  // const RankTwoTensor dxim1dE_tensor = dxidE_tensor * (-1.0 / (xi * xi));

  // // Compute terms for dSedE
  // const Real lambda_term = _lambda_const[_qp] - _damaged_modulus[_qp] / xi;
  // const Real shear_term = 2.0 * _shear_modulus[_qp] - _damaged_modulus[_qp] * xi;

  // RankFourTensor term_se1 = identity.outerProduct(-_damaged_modulus[_qp] * I1 * dxim1dE_tensor);
  // RankFourTensor term_se2 = identity.outerProduct(identity) * lambda_term;
  // RankFourTensor term_se3 = Ee.outerProduct(-_damaged_modulus[_qp] * dxidE_tensor);
  // RankFourTensor term_se4 = RankFourTensor(RankFourTensor::initIdentityFour) * shear_term;

  // RankFourTensor dSedE = term_se1 + term_se2 + term_se3 + term_se4;

  // // Compute terms for dSbdE
  // const Real coeff2_b = 2.0 * _a2[_qp] + _a1[_qp] / xi + 3.0 * _a3[_qp] * xi;
  // const Real coeff4_b = 2.0 * _a0[_qp] + _a1[_qp] * xi - _a3[_qp] * xi * xi * xi;

  // RankFourTensor term_b1 = identity.outerProduct((_a1[_qp] * dxim1dE_tensor + 3 * _a3[_qp] * dxidE_tensor) * I1);
  // RankFourTensor term_b2 = identity.outerProduct(identity) * coeff2_b;
  // RankFourTensor term_b3 = Ee.outerProduct(_a1[_qp] * dxidE_tensor - _a3[_qp] * 3 * xi * xi * dxidE_tensor);
  // RankFourTensor term_b4 = RankFourTensor(RankFourTensor::initIdentityFour) * coeff4_b;

  // RankFourTensor dSbdE = term_b1 + term_b2 + term_b3 + term_b4;

  // // Combine and assign tangent
  // tangent = (1.0 - _B_breakagevar[_qp]) * dSedE + _B_breakagevar[_qp] * dSbdE;  

}

// void
// ComputeLagrangianDamageBreakageStressPK2ModifiedFlowRule::computeQpPK2Stress()
// {
//   // Tolerances and maximum iterations for the predictor-corrector loop
//   const Real tol_theta = 1e-6;
//   const Real tol_tau   = 1e-6;
//   const Real tol_Dp    = 1e-6;
//   const unsigned int max_iter = 50;
  
//   // Old state (θ_old) from previous time step
//   const Real theta_old = _Theta_old[_qp];
//   // Initial guess for theta is the old value
//   Real theta_new = theta_old;
  
//   // Set initial guess for deviatoric stress from stored old value
//   RankTwoTensor tau_new = _Tau_old[_qp];

//   // Declare kinematics variables outside the loop so they remain in scope later
//   RankTwoTensor Fp_updated, Fe, Ee, Ep, E;
//   Real I1 = 0.0, I2 = 0.0, xi = 0.0;
//   RankTwoTensor Dp_new, Dp_old; // For comparing Dp updates

//   // Initialize Dp_old to zero for first iteration comparison.
//   Dp_old.zero();

//   // Loop to update theta and plastic strain rate concurrently with tau
//   for (unsigned int iter = 0; iter < max_iter; iter++)
//   {
//     // --- Compute Plastic Strain Rate Dp using current tau guess ---
//     Dp_new = _C_g[_qp] *
//              std::pow(_B_breakagevar_old[_qp], _m1[_qp]) *
//              tau_new *
//              std::pow(theta_new / _const_theta_o_mat[_qp], -1.0 * _const_B_mat[_qp] / _const_A_mat[_qp]);
    
//     // --- Compute Equivalent Plastic Strain Rate ---
//     Real Dp_eq = 0.0;
//     for (unsigned int i = 0; i < 3; i++){
//       for (unsigned int j = 0; j < 3; j++){
//         Dp_eq += Dp_new(i,j) * Dp_new(i,j);
//       }
//     }
//     Dp_eq = std::sqrt((2.0/3.0) * Dp_eq);
    
//     // --- Update Theta using Implicit Time Integration ---
//     Real dtheta_dt = 1.0 - Dp_eq * theta_new; 
//     Real theta_updated = theta_old + _dt * dtheta_dt;
    
//     // --- Compute Updated Kinematics & Stress using new Dp ---
//     // Update plastic deformation gradient
//     RankTwoTensor Cp = RankTwoTensor::Identity() - Dp_new * _dt;
//     Fp_updated = Cp.inverse() * _Fp_old[_qp];
    
//     // Compute elastic deformation gradient Fe
//     Fe = _F[_qp] * Fp_updated.inverse();
//     // Compute elastic strain: Ee = 0.5*(Fe^T * Fe - I)
//     Ee = 0.5 * (Fe.transpose() * Fe - RankTwoTensor::Identity());
    
//     // Compute invariants
//     I1 = Ee.trace();
//     I2 = 0.0;
//     for (unsigned int i = 0; i < 3; ++i){
//       for (unsigned int j = 0; j < 3; ++j){
//         I2 += Ee(i,j) * Ee(i,j);
//       }
//     }
    
//     // Compute xi (add a small number if necessary to avoid singularity)
//     xi = I1 / (std::sqrt(I2));
//     if (std::isnan(xi)) { xi = -std::sqrt(3); }
    
//     // --- Compute Stress Components ---
//     RankTwoTensor sigma_s = (_lambda_const[_qp] - _damaged_modulus[_qp] / xi) * I1 * RankTwoTensor::Identity() +
//                              (2 * _shear_modulus[_qp] - _damaged_modulus[_qp] * xi) * Ee;
//     RankTwoTensor sigma_b = (2 * _a2[_qp] + _a1[_qp] / xi + 3 * _a3[_qp] * xi) * I1 * RankTwoTensor::Identity() +
//                              (2 * _a0[_qp] + _a1[_qp] * xi - _a3[_qp] * std::pow(xi, 3)) * Ee;
//     RankTwoTensor sigma_total = (1 - _B_breakagevar[_qp]) * sigma_s + _B_breakagevar[_qp] * sigma_b;
    
//     // Save total stress (if needed)
//     _S[_qp] = sigma_total;
    
//     /* Compute plastic stress */
//     _Tp[_qp] = Fe.transpose() * Fe * sigma_total;
    
//     // --- Update Deviatoric Stress using updated plastic stress ---
//     RankTwoTensor tau_updated = _Tp[_qp] - 0.3333 * (_Tp[_qp].trace()) * RankTwoTensor::Identity();
    
//     // --- Convergence Check ---
//     // Check change in theta, tau, and Dp between iterations.
//     bool converged_theta = (std::fabs(theta_updated - theta_new) < tol_theta);
//     bool converged_tau   = ((tau_updated - tau_new).norm() < tol_tau);
//     bool converged_Dp    = ((Dp_new - Dp_old).norm() < tol_Dp);

//     // if (converged_theta && converged_tau && converged_Dp)
//     if (converged_theta)
//     {
//       theta_new = theta_updated;
//       _Theta[_qp] = theta_new;
//       _Dp[_qp] = Dp_new;
//       tau_new = tau_updated;
//       break;
//     }
    
//     // Update guesses for next iteration:
//     theta_new = theta_updated;
//     _Theta[_qp] = theta_new;
//     _Dp[_qp] = Dp_new;
//     tau_new = tau_updated;
//     Dp_old = Dp_new;
    
//     if (iter == max_iter - 1)
//     {
//       mooseWarning("Failed to converge in the predictor-corrector loop!");
//     }
//   }

//   // Save final updated tau state
//   _Tau[_qp] = tau_new;

//   /* Compute tangent */
//   RankFourTensor tangent;
//   computeQpTangentModulus(tangent, I1, I2, xi, Ee);

//   // Save tangent modulus
//   _C[_qp] = tangent;

//   /* Save other parameters */
//   _Fp[_qp] = Fp_updated;
//   _Fe[_qp] = Fe;
//   _Ee[_qp] = Ee;
//   // Make sure to compute E if needed before saving it:
//   // e.g., E = Fp_updated.transpose() * Ee * Fp_updated + Ep; // if Ep computed earlier
//   _E[_qp]  = E;
//   _I1[_qp] = I1;
//   _I2[_qp] = I2;
//   _xi[_qp] = xi;
// }