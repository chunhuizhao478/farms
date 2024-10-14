//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ComputeLagrangianDamageBreakageStressPK2.h"


registerMooseObject("farmsApp", ComputeLagrangianDamageBreakageStressPK2);

InputParameters
ComputeLagrangianDamageBreakageStressPK2::validParams()
{
  InputParameters params = ComputeLagrangianStressPK1::validParams();
  return params;
}

ComputeLagrangianDamageBreakageStressPK2::ComputeLagrangianDamageBreakageStressPK2(const InputParameters & parameters)
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
  _C(declareProperty<RankFourTensor>(_base_name + "pk2_jacobian")),
  _lambda_const(getMaterialProperty<Real>("lambda_const")),
  _shear_modulus(getMaterialProperty<Real>("shear_modulus")),
  _damaged_modulus(getMaterialProperty<Real>("damaged_modulus")),
  _B_breakagevar(getMaterialProperty<Real>("B_damagedvar")),
  _B_breakagevar_old(getMaterialPropertyOldByName<Real>("B_damagedvar")),
  _Tau_old(getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "deviatroic_stress")),
  _Fp_old(getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "plastic_deformation_gradient")),
  _F_old(getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "deformation_gradient")),
  _Ep_old(getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "plastic_strain")),
  _C_g(getMaterialProperty<Real>("C_g")),
  _m1(getMaterialProperty<Real>("m1")),
  _m2(getMaterialProperty<Real>("m2")),
  _a0(getMaterialProperty<Real>("a0")),
  _a1(getMaterialProperty<Real>("a1")),
  _a2(getMaterialProperty<Real>("a2")),
  _a3(getMaterialProperty<Real>("a3")),
  _dim(_mesh.dimension())
{
}

//Rules:See https://github.com/idaholab/moose/discussions/19450
//Only the object that declares the material property can assign values to it.
//Objects can request material properties, gaining read-only access to their values.
//When any object (including the object that declares it) requests the old value of a material property, that property becomes "stateful".
//All stateful material properties must be initialized within the initQpStatefulProperties call. 
//
void
ComputeLagrangianDamageBreakageStressPK2::initQpStatefulProperties()
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
}

void
ComputeLagrangianDamageBreakageStressPK2::computeQpPK1Stress()
{
  // PK2 update
  computeQpPK2Stress();

  // Compute Jp, Fp^{-1} //depends on the plastic deformation rate, here no volumetric strain upodate, Jp = 1 (checked)
  Real Jp = _Fp[_qp].det();
  //save
  _Jp[_qp] = Jp;

  RankTwoTensor Fpinv = _Fp[_qp].inverse();

  // Compute Fp_dot, F_dot
  // Here we approximate the rate by first-order, not sure if this is sufficient for varying time steps
  // currently MOOSE don't support getMaterialPropertyDot
  RankTwoTensor Fp_dot;
  RankTwoTensor F_dot;
  for (unsigned int i = 0; i < _dim; i++){
    for (unsigned int j = 0; j < _dim; j++){
      Fp_dot(i,j) += ( _Fp[_qp](i,j) - _Fp_old[_qp](i,j) ) / _dt;
      F_dot(i,j)  += (  _F[_qp](i,j) -  _F_old[_qp](i,j) ) / _dt;
    }
  }

  // Compute delta function
  auto delta = [](int i, int j) -> Real {
    return (i == j) ? 1.0 : 0.0;
  };  

  // Compute dFpdF //need to confirm
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

    //initialize value
    Real dFedF_val = delta(i,k) * Fpinv(l,m);

    //here apply summations to {h,r}
    for (unsigned int h = 0; h < _dim; h++){
      for (unsigned int r = 0; r < _dim; r++){
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
    for (unsigned int m = 0; m < _dim; m++){
      dEdF_val += 0.5 * ( dFedF(m,p,k,l) * _Fe[_qp](m,q) + _Fe[_qp](m,p) * dFedF(m,q,k,l) );
    }

    return dEdF_val;

  };

  // Compute dFpmdF
  auto dFpmdF = [&](int j, int n, int k, int l) -> Real {

    //initialize value
    Real dFpmdF_val = 0.0;

    //here apply summations to {i,a}
    for (unsigned int i = 0; i < _dim; i++){
      for (unsigned int m = 0; m < _dim; m++){
        dFpmdF_val += -1.0 * Fpinv(j,i) * dFpdF(i,m,k,l) * Fpinv(m,n);
      }
    }

    return dFpmdF_val;

  };

  //Compute pk_jacobian
  RankFourTensor pk_jacobian_val;
  pk_jacobian_val.zero();  // Make sure the tensor starts with zero values
  for (unsigned int i = 0; i < _dim; i++){
    for (unsigned int j = 0; j < _dim; j++){
      for (unsigned int k = 0; k < _dim; k++){
        for (unsigned int l = 0; l < _dim; l++){
          for (unsigned int m = 0; m < _dim; m++){
            // First term: dFedF(i,m,k,l) * S(m,n) * Fpinv(j,n)
            for (unsigned int n = 0; n < _dim; n++){
              pk_jacobian_val(i,j,k,l) += dFedF(i,m,k,l) * _S[_qp](m,n) * Fpinv(j,n);
            }
            
            // Second term: Fe(i,m) * C(m,n,p,q) * dEdF(p,q,k,l) * Fpinv(j,n)
            for (unsigned int n = 0; n < _dim; n++){
              for (unsigned int p = 0; p < _dim; p++){
                for (unsigned int q = 0; q < _dim; q++){
                  pk_jacobian_val(i,j,k,l) += _Fe[_qp](i,m) * _C[_qp](m,n,p,q) * dEdF(p,q,k,l) * Fpinv(j,n);
                }
              }
            }
            
            // Third term: Fe(i,m) * S(m,n) * dFpmdF(j,n,k,l)
            for (unsigned int n = 0; n < _dim; n++){
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
    _pk1_stress[_qp] = Jp * _Fe[_qp] * _S[_qp] * Fpinv.transpose();

    // Compute pk1 jacobian
    _pk1_jacobian[_qp] = pk_jacobian_val;

    //if there is no plastic 
    // _pk1_stress[_qp] = _F[_qp] * _S[_qp];
    // usingTensorIndices(i_, j_, k_, l_);
    // RankFourTensor dE =
    //     0.5 * (RankTwoTensor::Identity().times<i_, l_, j_, k_>(_F[_qp].transpose()) +
    //            _F[_qp].transpose().times<i_, k_, j_, l_>(RankTwoTensor::Identity()));

    // _pk1_jacobian[_qp] = RankTwoTensor::Identity().times<i_, k_, j_, l_>(_S[_qp].transpose()) +
    //                      (_C[_qp] * dE).singleProductI(_F[_qp]);

  }
  else
  {
    mooseError("Must selection 'large_kinematics' option!");
  }
}

void
ComputeLagrangianDamageBreakageStressPK2::computeQpPK2Stress()
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
  for (unsigned int i = 0; i < _dim; ++i){
    for (unsigned int j = 0; j < _dim; ++j){
      I2 += Ee(i,j) * Ee(i,j);
    }
  }

  /* Compute xi */
  //here we may need to add small number to avoid singularity
  Real xi = I1 / (std::sqrt(I2));

  /* Compute stress */
  RankTwoTensor sigma_s = (_lambda_const[_qp] - _damaged_modulus[_qp] / xi) * I1 * RankTwoTensor::Identity() + (2 * _shear_modulus[_qp] - _damaged_modulus[_qp] * xi) * Ee;
  RankTwoTensor sigma_b = (2 * _a2[_qp] + _a1[_qp] / xi + 3 * _a3[_qp] * xi) * I1 * RankTwoTensor::Identity() + (2 * _a0[_qp] + _a1[_qp] * xi - _a3[_qp] * std::pow(xi, 3)) * Ee;
  RankTwoTensor sigma_total = (1 - _B_breakagevar[_qp]) * sigma_s + _B_breakagevar[_qp] * sigma_b;

  //save
  _Ep[_qp] = Ep;
  _S[_qp] = sigma_total;

  //compute deviatroic stress tensor //save
  _Tau[_qp] = sigma_total - 1.0 / 3.0 * ( sigma_total.trace() ) * RankTwoTensor::Identity();

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

RankTwoTensor
ComputeLagrangianDamageBreakageStressPK2::computeQpFp()
{
  //Apply power operation on every element of Tau
  //let's assume m2 = 1, and not apply pow on its elements
  RankTwoTensor Tau_old_power_m2 = _Tau_old[_qp];

  //Compute Plastic Deformation Rate Tensor Dp at t_{n+1} using quantities from t_{n}
  RankTwoTensor Dp = _C_g[_qp] * std::pow(_B_breakagevar_old[_qp],_m1[_qp]) * Tau_old_power_m2; 

  // //Update Plastic Strain
  // _Ep[_qp] = _Ep_old[_qp] + Dp * _dt;

  //Compute Cp = I - Dp dt
  RankTwoTensor Cp = RankTwoTensor::Identity() - Dp * _dt;

  //Use Implicit Euler Integration, Update Fp
  RankTwoTensor Fp_updated = Cp.inverse() * _Fp_old[_qp];

  return Fp_updated;
}

void
ComputeLagrangianDamageBreakageStressPK2::computeQpTangentModulus(RankFourTensor & tangent, 
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
    return delta(k,l) * pow(I2,-0.5) - 0.5 * pow(I2,-1.5) * dI2dE(k,l) * I1;
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
  for (unsigned int i = 0; i < _dim; i++){
    for (unsigned int j = 0; j < _dim; j++){
      for (unsigned int k = 0; k < _dim; k++){
        for (unsigned int l = 0; l < _dim; l++){
          tangent(i,j,k,l) += dSdE(i,j,k,l);
        }
      }
    }
  }

}