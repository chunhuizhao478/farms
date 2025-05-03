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
  return params;
}

ComputeLagrangianDamageBreakageStressPK2Diffused::ComputeLagrangianDamageBreakageStressPK2Diffused(const InputParameters & parameters)
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
  //---------------------------------------------------------------------------------------------//
  // Add option to add dilatancy/compaction effect //Follow paper Section 7.1
  _eta(declareProperty<Real>(_base_name + "plastic_volume_change")),
  _dilatancy_function_beta(declareProperty<Real>(_base_name + "dilatancy_function_beta")),
  _shear_rate_nu(declareProperty<Real>(_base_name + "shear_rate_nu")),
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
  // Add option to add dilatancy/compaction effect //Follow paper Section 7.1
  _eta_old(getMaterialPropertyOldByName<Real>(_base_name + "plastic_volume_change")),
  //---------------------------------------------------------------------------------------------//
  _C_g(getMaterialProperty<Real>("C_g")),
  _m1(getMaterialProperty<Real>("m1")),
  _m2(getMaterialProperty<Real>("m2")),
  _a0(getMaterialProperty<Real>("a0")),
  _a1(getMaterialProperty<Real>("a1")),
  _a2(getMaterialProperty<Real>("a2")),
  _a3(getMaterialProperty<Real>("a3")),
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
  _shear_stress_perturbation(getMaterialProperty<Real>("shear_stress_perturbation"))
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
}

void
ComputeLagrangianDamageBreakageStressPK2Diffused::computeQpPK1Stress()
{
  // // PK2 update
  // computeQpPK2Stress();

  // // Compute Jp, Fp^{-1} //depends on the plastic deformation rate, here no volumetric strain upodate, Jp = 1 (checked)
  // _Jp[_qp] = _Fp[_qp].det();
  // RankTwoTensor Fpinv = _Fp[_qp].inverse();

  // //Compute deformation rate D
  // _D[_qp] = 0.5 * ( _velgrad_L[_qp] + _velgrad_L[_qp].transpose() );

  // // Compute Fp_dot, F_dot
  // // Here we approximate the rate by first-order, not sure if this is sufficient for varying time steps
  // // currently MOOSE don't support getMaterialPropertyDot
  // RankTwoTensor Fp_dot;
  // RankTwoTensor F_dot;
  
  // if (_use_vels_build_L_mat[_qp]){ //use true deformation rate

  //   for (unsigned int i = 0; i < 3; i++){
  //     for (unsigned int j = 0; j < 3; j++){
  //       for (unsigned int m = 0; m < 3; m++){
  //         Fp_dot(i,j) += _Dp[_qp](i,m) * _Fp[_qp](m,j);
  //         F_dot(i,j) += _D[_qp](i,m) * _F[_qp](m,j); 
  //       }
  //     }
  //   } 

  // }
  // else{ //finite difference approximation

  //   for (unsigned int i = 0; i < 3; i++){
  //     for (unsigned int j = 0; j < 3; j++){
  //         F_dot(i,j)  = (_F[_qp](i,j) - _F_old[_qp](i,j) ) / _dt; 
  //         //F_dot(i,j)  = (_F[_qp](i,j) - _F_old[_qp](i,j) ); 
  //         for (unsigned int m = 0; m < 3; m++){
  //           Fp_dot(i,j) += _Dp[_qp](i,m) * _Fp[_qp](m,j);
  //         }
  //     }
  //   }
  
  // }
 
  // // Compute delta function
  // auto delta = [](int i, int j) -> Real {
  //   return (i == j) ? 1.0 : 0.0;
  // };  

  // //Compute dFpdF //need to confirm
  // auto dFpdF = [&](int i, int j, int k, int l) -> Real {
  //   if (_dt == 0.0){ //Steady, set zero
  //     return 0.0;
  //   } 
  //   else{ //Transient
  //     if (Fp_dot(i,j) == 0.0 || F_dot(k,l) == 0.0 ){ //no change of viscoelastic dg, set zero
  //       return 0.0;
  //     }
  //     else{
  //       return Fp_dot(i,j)/F_dot(k,l);
  //     }
  //   }
  // };

  // // Compute dFedF
  // auto dFedF = [&](int i, int m, int k, int l) -> Real {

  //   // Since dFp/dF is zero (Fp is explicit), dFedF simplifies
  //   Real dFedF_val = delta(i,k) * Fpinv(l,m);

  //   //here apply summations to {h,r}
  //   for (unsigned int h = 0; h < 3; h++){
  //     for (unsigned int r = 0; r < 3; r++){
  //       dFedF_val -= _Fe[_qp](i,h) * dFpdF(h,r,k,l) * Fpinv(r,m);
  //     }
  //   }    

  //   return dFedF_val;

  // }; 

  // // Compute dEdF
  // auto dEdF = [&](int p, int q, int k, int l) -> Real {
    
  //   //initialize value
  //   Real dEdF_val = 0.0;

  //   //here apply summations to {m}
  //   for (unsigned int m = 0; m < 3; m++){
  //     dEdF_val += 0.5 * ( dFedF(m,p,k,l) * _Fe[_qp](m,q) + _Fe[_qp](m,p) * dFedF(m,q,k,l) );
  //   }

  //   return dEdF_val;

  // };

  // // Compute dFpmdF
  // auto dFpmdF = [&](int j, int n, int k, int l) -> Real {

  //   //initialize value
  //   Real dFpmdF_val = 0.0;

  //   //here apply summations to {i,a}
  //   for (unsigned int i = 0; i < 3; i++){
  //     for (unsigned int m = 0; m < 3; m++){
  //       dFpmdF_val += -1.0 * Fpinv(j,i) * dFpdF(i,m,k,l) * Fpinv(m,n);
  //     }
  //   }

  //   return dFpmdF_val;

  // };

  // //Compute pk_jacobian
  // RankFourTensor pk_jacobian_val;
  // pk_jacobian_val.zero();  // Make sure the tensor starts with zero values
  // for (unsigned int i = 0; i < 3; i++){
  //   for (unsigned int j = 0; j < 3; j++){
  //     for (unsigned int k = 0; k < 3; k++){
  //       for (unsigned int l = 0; l < 3; l++){
  //         for (unsigned int m = 0; m < 3; m++)
  //         {
  //           // First term: dFedF(i,m,k,l) * S(m,n) * Fpinv(j,n)
  //           for (unsigned int n = 0; n < 3; n++){
  //             pk_jacobian_val(i,j,k,l) += dFedF(i,m,k,l) * _S[_qp](m,n) * Fpinv(j,n);
  //           }
            
  //           // Second term: Fe(i,m) * C(m,n,p,q) * dEdF(p,q,k,l) * Fpinv(j,n)
  //           for (unsigned int n = 0; n < 3; n++){
  //             for (unsigned int p = 0; p < 3; p++){
  //               for (unsigned int q = 0; q < 3; q++){
  //                 pk_jacobian_val(i,j,k,l) += _Fe[_qp](i,m) * _C[_qp](m,n,p,q) * dEdF(p,q,k,l) * Fpinv(j,n);
  //               }
  //             }
  //           }
            
  //           // Third term: Fe(i,m) * S(m,n) * dFpmdF(j,n,k,l)
  //           for (unsigned int n = 0; n < 3; n++){
  //             pk_jacobian_val(i,j,k,l) += _Fe[_qp](i,m) * _S[_qp](m,n) * dFpmdF(j,n,k,l);
  //           }
  //         }
  //       }
  //     }
  //   }
  // }

  // // Complicated wrapping from PK2 to PK1, see documentation on overleaf
  // if (_large_kinematics)
  // {
  //   //if there is plastic
  //   // Compute pk1 stress
  //   _pk1_stress[_qp] = _Fe[_qp] * _S[_qp] * Fpinv.transpose();

  //   // Compute pk1 jacobian
  //   _pk1_jacobian[_qp] = pk_jacobian_val;

  //   //if there is no plastic 
  //   // _pk1_stress[_qp] = _F[_qp] * _S[_qp];
  //   // usingTensorIndices(i_, j_, k_, l_);
  //   // RankFourTensor dE =
  //   //     0.5 * (RankTwoTensor::Identity().times<i_, l_, j_, k_>(_F[_qp].transpose()) +
  //   //            _F[_qp].transpose().times<i_, k_, j_, l_>(RankTwoTensor::Identity()));

  //   // _pk1_jacobian[_qp] = RankTwoTensor::Identity().times<i_, k_, j_, l_>(_S[_qp].transpose()) +
  //   //                      (_C[_qp] * dE).singleProductI(_F[_qp]);

  // }
  // else
  // {
  //   mooseError("Must selection 'large_kinematics' option!");
  // }

  //--------------------------------------------------------------------------
  // PK2 update
  computeQpPK2Stress();

  // // Compute Jp and the inverse of Fp
  // if (_add_dilatancy_compaction_anand_mat[_qp]){
  //   //here we assume exponential dependence of the plastic volume change eta: 
  //   //eta = ln(J^p), J^p = exp(eta)
  //   _Jp[_qp] = std::exp(_eta[_qp]);
  // }
  // else{
  //   _Jp[_qp] = _Fp[_qp].det();
  // }

  _Jp[_qp] = _Fp[_qp].det();
  
  RankTwoTensor Fpinv = _Fp[_qp].inverse();

  //Compute deformation rate D
  // _D[_qp] = 0.5 * ( _velgrad_L[_qp] + _velgrad_L[_qp].transpose() );

  // Compute Fp_dot, F_dot
  // Here we approximate the rate by first-order, not sure if this is sufficient for varying time steps
  // currently MOOSE don't support getMaterialPropertyDot
  RankTwoTensor Fp_dot; Fp_dot.zero();
  RankTwoTensor F_dot; F_dot.zero();
  
  // if (_use_vels_build_L_mat[_qp]){ //use true deformation rate

  //   for (unsigned int i = 0; i < 3; i++){
  //     for (unsigned int j = 0; j < 3; j++){
  //       for (unsigned int m = 0; m < 3; m++){
  //         Fp_dot(i,j) += _Dp[_qp](i,m) * _Fp[_qp](m,j);
  //         F_dot(i,j) += _D[_qp](i,m) * _F[_qp](m,j); 
  //       }
  //     }
  //   } 

  // }
  // else{ //finite difference approximation

  //   for (unsigned int i = 0; i < 3; i++){
  //     for (unsigned int j = 0; j < 3; j++){
  //         //F_dot(i,j)  = (_F[_qp](i,j) - _F_old[_qp](i,j) ) / _dt; 
  //         F_dot(i,j)  = (_F[_qp](i,j) - _F_old[_qp](i,j) ); 
  //         for (unsigned int m = 0; m < 3; m++){
  //           Fp_dot(i,j) += _Dp[_qp](i,m) * _Fp[_qp](m,j);
  //         }
  //     }
  //   }
  
  // }

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
  // PK2-to-PK1 wrapping: using large kinematics formulation.
  if (_large_kinematics)
  {
    // Compute PK1 stress: P = Jp * Fe * S * (Fpinv)^T
    _pk1_stress[_qp] = _Jp[_qp] * _Fe[_qp] * _S[_qp] * Fpinv.transpose();
    // Here we assume Jp = 1
    //_pk1_stress[_qp] = _Fe[_qp] * _S[_qp] * Fpinv.transpose();
    // Assign the computed consistent tangent operator
    _pk1_jacobian[_qp] = pk_jacobian_val;
  }
  else
  {
    mooseError("Must select 'large_kinematics' option!");
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
  //RankTwoTensor E = 0.5 * (_F[_qp].transpose() * _F[_qp] - RankTwoTensor::Identity());

  // //----------------------------------------------------------------------------//
  // // Add check to ensure Ep is always not exceeding E
  // Real Eij = 0.0;
  // Real Epij = 0.0;
  // for (unsigned int i = 0; i < 3; ++i){
  //   for (unsigned int j = 0; j < 3; ++j){
  //     Eij  = E(i,j);
  //     Epij = Ep(i,j);
  //     if (Eij > 0.0 && Epij > Eij){
  //       Epij = Eij;
  //     }
  //     else if (Eij < 0.0 && Epij < Eij){
  //       Epij = Eij;
  //     }
  //     else if (Eij == 0.0 && Epij != 0.0){
  //       Epij = 0.0;
  //     }
  //     Ep(i,j) = Epij;
  //   }
  // }

  // // Update Ee
  // // 1) Precompute inverses of Fp
  // RankTwoTensor Fp_inv  = Fp_updated.inverse();
  // RankTwoTensor Fp_invT = Fp_inv.transpose();

  // // 2) Compute Ee_corrected = Fp_invT * (E - Ep) * Fp_inv
  // RankTwoTensor E_minus_Ep = E - Ep;
  // Ee = Fp_invT * E_minus_Ep * Fp_inv;
  // Ee = 0.5 * (Ee + Ee.transpose());
  // //----------------------------------------------------------------------------//

  //Compute shear perturbation
  //----------------------------------------------------------------------------//
  Real shear_strain_perturbation = 0.0;
  if (_shear_stress_perturbation[_qp] != 0){
    shear_strain_perturbation = _shear_stress_perturbation[_qp] / (2 * _shear_modulus[_qp]);
    Ee(0,1) += shear_strain_perturbation;
    Ee(1,0) += shear_strain_perturbation;
  }
  //----------------------------------------------------------------------------//

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

  // add structral stress
  for (unsigned int i = 0; i < 3; ++i){
    for (unsigned int j = 0; j < 3; ++j){
      sigma_total(i,j) -= _structural_stress_coefficient[_qp] * _grad_alpha_damagedvar[_qp](i) * _grad_alpha_damagedvar[_qp](j);
    }
  }

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

  //compute deviatoric strain rate
  computeDeviatroicStrainRateTensor();

}

RankTwoTensor
ComputeLagrangianDamageBreakageStressPK2Diffused::computeQpFp()
{
  //NOT FINISHED
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

  //Get deviatroic stress direction
  RankTwoTensor N; N.zero();
  // Epsilon to avoid division by zero
  if (Tau_eq != 0.0){
    //Compute deviatroic stress direction
    for (unsigned int p = 0; p < 3; p++){
      for (unsigned int q = 0; q < 3; q++){
        N(p,q) = Tau_old(p,q) / Tau_eq;
      }
    }
  }

  //Define equvialent shear rate nu
  if (_use_state_var_evolution_mat[_qp]){
    computestatedependentDp();
  }
  else{
    _shear_rate_nu[_qp] = _C_g[_qp] * std::pow(_B_breakagevar_old[_qp],_m1[_qp]) * std::pow(Tau_eq,_m2[_qp]);
  }

  //Compute Plastic Deformation Rate Tensor Dp at t_{n+1} using quantities from t_{n}
  RankTwoTensor Dp = _shear_rate_nu[_qp] * N; 

  // //Add dilatancy/compaction effect
  // if (_add_dilatancy_compaction_anand_mat[_qp]){
    
  //   //Compute dilatancy function beta //_dilatancy_function_beta[_qp]
  //   computedilatancyfunction();

  //   //Update Plastic Deformation Rate Tensor
  //   Dp = Dp + _dilatancy_function_beta[_qp] * _shear_rate_nu[_qp] * RankTwoTensor::Identity();

  //   //Update plastic volume change
  //   computeplasticvolumechange();

  // }

  // //Update Plastic Strain
  // _Ep[_qp] = _Ep_old[_qp] + Dp * _dt;

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

  // //define functions for derivatives

  // //delta function
  // auto delta = [](int i, int j) -> Real {
  //   return (i == j) ? 1.0 : 0.0;
  // };

  // //dI1_dE_{kl}
  // auto dI1dE = [&](int k, int l) -> Real {
  //   return delta(k,l);
  // };

  // //dI2_dE_{kl}
  // auto dI2dE = [&](int k, int l) -> Real {
  //   return 2 * Ee(k,l);
  // };

  // //dxi_dE_{kl}
  // auto dxidE = [&](int k, int l) -> Real {
    
  //   // Epsilon to avoid division by zero
  //   const Real epsilon = 1e-12;
  //   // Adjust I2 if necessary
  //   Real adjusted_I2 = I2;
  //   if (I2 <= epsilon) {
  //     //mooseWarning("I2 is zero or too small (I2 = ", I2, "), adjusting to epsilon.");
  //     adjusted_I2 = epsilon;
  //   }

  //   Real dxidE = 0.5 * pow(adjusted_I2,-1.5) * dI2dE(k,l) * I1;
  //   //mooseInfo("I1 = ", I1, ", I2 = ", I2);
  //   if (std::isnan(dxidE)){mooseError("dxidE");}
  //   return delta(k,l) * pow(adjusted_I2,-0.5) - 0.5 * pow(adjusted_I2,-1.5) * dI2dE(k,l) * I1;
  // };

  // //dE_{ij}_dE_{kl}
  // auto dEdE = [&](int i, int j, int k, int l) -> Real {
  //   return delta(i,k) * delta(j,l);
  //   //return 0.5 * ( delta(i,k) * delta(j,l) + delta(i,l) * delta(j,k) ); //its symmetric form
  // };

  // //dxi^{-1}_dE_{kl}
  // auto dxim1dE = [&](int k, int l) -> Real {
  //   return -1.0 * pow(xi,-2.0) * dxidE(k,l);
  // };

  // //dxi^3_dE_{kl}
  // auto dxi3dE = [&](int k, int l) -> Real {
  //   return 3 * pow(xi,2) * dxidE(k,l);
  // };

  // //dSe_{ij}_dE_{kl}
  // auto dSedE = [&](int i, int j, int k, int l) -> Real {
  //   Real dSedE_components = (- _damaged_modulus[_qp] * dxim1dE(k,l) ) * I1 * delta(i,j);
  //   dSedE_components += ( _lambda_const[_qp] - _damaged_modulus[_qp] / xi ) * dI1dE(k,l) * delta(i,j);
  //   dSedE_components += (- _damaged_modulus[_qp] * dxidE(k,l) ) * Ee(i,j);
  //   dSedE_components += ( 2 * _shear_modulus[_qp] - _damaged_modulus[_qp] * xi ) * dEdE(i,j,k,l);
  //   if (std::isnan(dSedE_components)){mooseError("dSedE_components");}
  //   return dSedE_components;
  // };

  // //dSb_{ij}_dE_{kl}
  // auto dSbdE = [&](int i, int j, int k, int l) -> Real {
  //   Real dSbdE_components = ( _a1[_qp] * dxim1dE(k,l) + 3 * _a3[_qp] * dxidE(k,l) ) * I1 * delta(i,j);
  //   dSbdE_components += ( 2 * _a2[_qp] + _a1[_qp] / xi + 3 * _a3[_qp] * xi ) * dI1dE(k,l) * delta(i,j);
  //   dSbdE_components += ( _a1[_qp] * dxidE(k,l) - _a3[_qp] * dxi3dE(k,l) ) * Ee(i,j);
  //   dSbdE_components += ( 2 * _a0[_qp] + _a1[_qp] * xi - _a3[_qp] * pow(xi,3) ) * dEdE(i,j,k,l);
  //   return dSbdE_components;
  // };

  // //dS_{ij}_dE_{kl}
  // auto dSdE = [&](int i, int j, int k, int l) -> Real {
  //   return (1 - _B_breakagevar[_qp]) * dSedE(i,j,k,l) + _B_breakagevar[_qp] * dSbdE(i,j,k,l);
  // };

  // // Compute tangent modulus C
  // for (unsigned int i = 0; i < 3; i++){
  //   for (unsigned int j = 0; j < 3; j++){
  //     for (unsigned int k = 0; k < 3; k++){
  //       for (unsigned int l = 0; l < 3; l++){
  //         if (std::isnan(dSdE(i,j,k,l))){mooseError("encounter nan error: dSdE(i,j,k,l)");}
  //         tangent(i,j,k,l) += dSdE(i,j,k,l);
  //       }
  //     }
  //   }
  // }

  const Real adjusted_I2 = (I2 <= 1e-12) ? 1e-12 : I2;
  const RankTwoTensor identity = RankTwoTensor::Identity();

  // Precompute dxidE tensor
  RankTwoTensor dxidE_tensor;
  for (unsigned int k = 0; k < 3; ++k)
    for (unsigned int l = 0; l < 3; ++l)
      dxidE_tensor(k, l) = (identity(k, l) * adjusted_I2 - I1 * Ee(k, l)) / std::pow(adjusted_I2, 1.5);

  const RankTwoTensor dxim1dE_tensor = dxidE_tensor * (-1.0 / (xi * xi));

  // Compute terms for dSedE
  const Real lambda_term = _lambda_const[_qp] - _damaged_modulus[_qp] / xi;
  const Real shear_term = 2.0 * _shear_modulus[_qp] - _damaged_modulus[_qp] * xi;

  RankFourTensor term_se1 = identity.outerProduct(-_damaged_modulus[_qp] * I1 * dxim1dE_tensor);
  RankFourTensor term_se2 = identity.outerProduct(identity) * lambda_term;
  RankFourTensor term_se3 = Ee.outerProduct(-_damaged_modulus[_qp] * dxidE_tensor);
  RankFourTensor term_se4 = RankFourTensor(RankFourTensor::initIdentityFour) * shear_term;

  RankFourTensor dSedE = term_se1 + term_se2 + term_se3 + term_se4;

  // Compute terms for dSbdE
  const Real coeff2_b = 2.0 * _a2[_qp] + _a1[_qp] / xi + 3.0 * _a3[_qp] * xi;
  const Real coeff4_b = 2.0 * _a0[_qp] + _a1[_qp] * xi - _a3[_qp] * xi * xi * xi;

  RankFourTensor term_b1 = identity.outerProduct((_a1[_qp] * dxim1dE_tensor + 3 * _a3[_qp] * dxidE_tensor) * I1);
  RankFourTensor term_b2 = identity.outerProduct(identity) * coeff2_b;
  RankFourTensor term_b3 = Ee.outerProduct(_a1[_qp] * dxidE_tensor - _a3[_qp] * 3 * xi * xi * dxidE_tensor);
  RankFourTensor term_b4 = RankFourTensor(RankFourTensor::initIdentityFour) * coeff4_b;

  RankFourTensor dSbdE = term_b1 + term_b2 + term_b3 + term_b4;

  // Combine and assign tangent
  tangent = (1.0 - _B_breakagevar[_qp]) * dSedE + _B_breakagevar[_qp] * dSbdE;  

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

// //Compute dilatancy function beta
// void
// ComputeLagrangianDamageBreakageStressPK2Diffused::computedilatancyfunction()
// {
//   _dilatancy_function_beta[_qp] = _anand_param_go_mat[_qp] * std::pow( 1 - _eta[_qp] / _anand_param_eta_cv_mat[_qp] , _anand_param_p_mat[_qp] );
// }

// //Compute plastic volume change eta
// void
// ComputeLagrangianDamageBreakageStressPK2Diffused::computeplasticvolumechange()
// {
//   _eta[_qp] = _eta_old[_qp] + _dilatancy_function_beta[_qp] * _shear_rate_nu[_qp] * _dt;
// }