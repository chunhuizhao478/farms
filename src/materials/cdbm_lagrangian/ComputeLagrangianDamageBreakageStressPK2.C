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
  _Fe(declareProperty<RankTwoTensor>(_base_name + "elastic_deformation_gradient")),
  _Tau(declareProperty<RankTwoTensor>(_base_name + "deviatroic_stress")),
  _Ee(declareProperty<RankTwoTensor>(_base_name + "green_lagrange_elastic_strain")),
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
  _I1[_qp] = 0.0;
  _I2[_qp] = 0.0;
  _xi[_qp] = -sqrt(3);
  _S[_qp].zero();
  _C[_qp].zero();
}

void
ComputeLagrangianDamageBreakageStressPK2::computeQpPK1Stress()
{
  
  RankTwoTensor zeroranktwo(0,0,0,0,0,0,0,0,0);
  _S[_qp] = zeroranktwo;
  _C[_qp].zero();

  // PK2 update
  computeQpPK2Stress();

  // Compute Jp, Fp^{-1}
  Real Jp = _Fp[_qp].det();
  RankTwoTensor Fpinv = _Fp[_qp].inverse();

  // Compute Fp_dot, F_dot
  // Here we approximate the rate by first-order, not sure if this is sufficient for varying time steps
  // currently MOOSE don't support getMaterialPropertyDot
  RankTwoTensor Fp_dot = ( _Fp[_qp] - _Fp_old[_qp] ) / _dt;
  RankTwoTensor F_dot  = (  _F[_qp] -  _F_old[_qp] ) / _dt;

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
      if (Fp_dot(i,j) == 0.0){ //no change of viscoelastic dg, set zero
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
    Real dFedF_val = 0.0;

    //here apply summations to {j}
    for (unsigned int j = 0; j < _dim; j++){
      dFedF_val += ( delta(i,k) * delta(j,l) - _Fe[_qp](i,m) * dFpdF(m,j,k,l) ) * Fpinv(j,m);
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
  auto dFpmdF = [&](int j, int m, int k, int l) -> Real {

    //initialize value
    Real dFpmdF_val = 0.0;

    //here apply summations to {i}
    for (unsigned int i = 0; i < _dim; i++){
      dFpmdF_val += -1 * Fpinv(i,m) * dFpdF(m,i,k,l) * Fpinv(j,m);
    }

    return dFpmdF_val;

  };

  // Compute dPdF
  auto dPdF = [&](int i, int j, int k, int l) -> Real {

    //initialize value
    Real dPdF_val = 0.0;

    // std::cout << "S in dPdF: " << std::endl;
    // S.print();

    // std::cout << "C in dPdF: " << std::endl;
    // C.print();

    for (unsigned int m = 0; m < _dim; m++) {
      for (unsigned int n = 0; n < _dim; n++) {
        for (unsigned int p = 0; p < _dim; p++) {
          for (unsigned int q = 0; q < _dim; q++) {
            // Output values being used in the computation //dFedF, dEdF, dFpmdF
            // std::cout << "dFedF(" << i << "," << m << "," << k << "," << l << "): " << dFedF(i, m, k, l) << std::endl;
            // std::cout << "_S[_qp](" << m << "," << n << "): " << _S[_qp](m, n) << std::endl;
            // std::cout << "Fpinv(" << j << "," << n << "): " << Fpinv(j, n) << std::endl;
            // std::cout << "_Fe[_qp](" << i << "," << m << "): " << _Fe[_qp](i, m) << std::endl;
            // std::cout << "_C[_qp](" << m << "," << n << "," << p << "," << q << "): " << _C[_qp](m, n, p, q) << std::endl;
            // std::cout << "dEdF(" << p << "," << q << "," << k << "," << l << "): " << dEdF(p, q, k, l) << std::endl;
            // std::cout << "dFpmdF(" << j << "," << n << "," << k << "," << l << "): " << dFpmdF(j, n, k, l) << std::endl;
            
            // Check for NaN before adding to dPdF_val
            Real current_value = dFedF(i, m, k, l) * _S[_qp](m, n) * Fpinv(j, n) +
                                _Fe[_qp](i, m) * _C[_qp](m, n, p, q) * dEdF(p, q, k, l) * Fpinv(j, n) +
                                _Fe[_qp](i, m) * _S[_qp](m, n) * dFpmdF(j, n, k, l);
            
            // if (std::isnan(current_value)) {
            //     std::cout << "NaN detected at indices: i=" << i << ", j=" << j << ", k=" << k << ", l=" << l << ", m=" << m << ", n=" << n << ", p=" << p << ", q=" << q << std::endl;
            // }
            
            dPdF_val += current_value;
          }
        }
      }
    }

    return dPdF_val;

  }; 

  // Compute pk_jacobian
  RankFourTensor pk_jacobian_val;
  for (unsigned int i = 0; i < _dim; i++){
    for (unsigned int j = 0; j < _dim; j++){
      for (unsigned int k = 0; k < _dim; k++){
        for (unsigned int l = 0; l < _dim; l++){
          pk_jacobian_val(i,j,k,l) = dPdF(i,j,k,l);
        }
      }
    }
  }

  // //check err //pk_jacobian_val has issue
  // std::cout<<"----------PK1-----------"<<std::endl;
  // std::cout<<"_Fe[_qp]"<<std::endl;
  // RankTwoTensor Fe_check = _Fe[_qp];
  // Fe_check.print();
  // std::cout<<"_S[_qp]"<<std::endl;
  // RankTwoTensor S_check = _S[_qp];
  // _S[_qp].print();
  // std::cout<<"Fpinv.transpose()"<<std::endl;
  // Fpinv.transpose().print();  
  // std::cout<<"pk_jacobian_val"<<std::endl;
  // pk_jacobian_val.print();
  // std::cout<<"-----------------------"<<std::endl;

  // Complicated wrapping from PK2 to PK1, see documentation on overleaf
  if (_large_kinematics)
  {
    // Compute pk1 stress
    _pk1_stress[_qp] = Jp * _Fe[_qp] * _S[_qp] * Fpinv.transpose();

    // Compute pk1 jacobian
    _pk1_jacobian[_qp] = pk_jacobian_val;

  }
  // Small deformations all are equivalent
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
  _S[_qp] = sigma_total;

  //compute deviatroic stress tensor //save
  _Tau[_qp] = sigma_total - 1.0 / 3.0 * ( sigma_total(0,0) + sigma_total(1,1) + sigma_total(2,2) ) * RankTwoTensor::Identity();

  /* Compute tangent */
  std::vector<Real> components(21);

  computeQpTangentModulus(components,I1,I2,xi,Ee);

  RankFourTensor tangent(components, RankFourTensor::symmetric21);

  //save
  _C[_qp] = tangent;

  /* Save other parameters */
  _Fp[_qp] = Fp_updated;
  _Fe[_qp] = Fe;
  _Ee[_qp] = Ee;
  _I1[_qp] = I1;
  _I2[_qp] = I2;
  _xi[_qp] = xi;

  // // //check err
  // std::cout<<"----------PK2-----------"<<std::endl;
  // for (unsigned int i = 0; i < _dim; ++i){
  //   for (unsigned int j = 0; j < _dim; ++j){
  //     if (std::isnan(sigma_s(i,j))){
  //       std::cout<<"find nan in sigma_s"<<std::endl;
  //     }
  //     if (std::isnan(sigma_b(i,j))){
  //       std::cout<<"find nan in sigma_b"<<std::endl;
  //     }
  //   }
  // }
  // if (std::isnan(_lambda_const[_qp])){ std::cout<<"find nan in _lambda_const[_qp]"<<std::endl; }else{ std::cout<<"_lambda_const[_qp]:"<<_lambda_const[_qp]<<std::endl;}
  // if (std::isnan(_shear_modulus[_qp])){ std::cout<<"find nan in _shear_modulus[_qp]"<<std::endl; }else{ std::cout<<"_shear_modulus[_qp]:"<<_shear_modulus[_qp]<<std::endl;}
  // if (std::isnan(_damaged_modulus[_qp])){ std::cout<<"find nan in _damaged_modulus[_qp]"<<std::endl; }else{ std::cout<<"_damaged_modulus[_qp]:"<<_damaged_modulus[_qp]<<std::endl;}
  // if (std::isnan(xi)){ std::cout<<"find nan in xi"<<std::endl; }else{ std::cout<<"xi:"<<xi<<std::endl;}
  // std::cout<<"Fp_updated"<<std::endl;
  // Fp_updated.print();
  // std::cout<<"Fe"<<std::endl;
  // Fe.print();
  // std::cout<<"Ee"<<std::endl;
  // Ee.print();  
  // std::cout<<"sigma_total"<<std::endl;
  // sigma_total.print();
  // std::cout<<"tangent"<<std::endl;
  // tangent.print();
  // std::cout<<"I1: "<<I1<<std::endl;
  // std::cout<<"I2: "<<I2<<std::endl;
  // std::cout<<"xi: "<<xi<<std::endl;
  // std::cout<<"-----------------------"<<std::endl;

}

RankTwoTensor
ComputeLagrangianDamageBreakageStressPK2::computeQpFp()
{
  //Apply power operation on every element of Tau
  RankTwoTensor Tau_old_power_m2;
  for (unsigned int i = 0; i < _dim; i++){
    for (unsigned int j = 0; j < _dim; j++){
      Tau_old_power_m2(i,j) = std::pow(_Tau_old[_qp](i,j), _m2[_qp]);
    }
  }

  //Compute Plastic Deformation Rate Tensor Dp at t_{n+1} using quantities from t_{n}
  RankTwoTensor Dp = _C_g[_qp] * std::pow(_B_breakagevar_old[_qp], _m1[_qp]) * Tau_old_power_m2; 

  //Compute Cp = I - Dp dt
  RankTwoTensor Cp = RankTwoTensor::Identity() - Dp * _dt;

  //Use Implicit Euler Integration, Update Fp
  RankTwoTensor Fp_updated = Cp.inverse() * _Fp_old[_qp];

  return Fp_updated;
}

void
ComputeLagrangianDamageBreakageStressPK2::computeQpTangentModulus(std::vector<Real>& tangent, 
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
    return 0.5 * ( delta(i,k) * delta(j,l) + delta(i,l) * delta(j,k) );
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
    Real dSedE_components = ( _lambda_const[_qp] - _damaged_modulus[_qp] * dxim1dE(k,l) ) * I1 * delta(i,j);
    dSedE_components += ( _lambda_const[_qp] - _damaged_modulus[_qp] / xi ) * dI1dE(k,l) * delta(i,j);
    dSedE_components += ( 2 * _shear_modulus[_qp] - _damaged_modulus[_qp] * dxidE(k,l) ) * Ee(i,j);
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
  //here we transform index 1,2,3 -> 0,1,2
  auto dSdE = [&](int i, int j, int k, int l) -> Real {
    int ind_i = i - 1; int ind_j = j - 1; int ind_k = k - 1; int ind_l = l - 1;
    return (1 - _B_breakagevar[_qp]) * dSedE(ind_i,ind_j,ind_k,ind_l) + _B_breakagevar[_qp] * dSbdE(ind_i,ind_j,ind_k,ind_l);
  };

  //fill in tangent vector
  //symmetric21 : C_ijkl = '1111 1122 1133 1123 1113 1112 2222 2233 2223 2213 2212 3333 3323 3313 3312 2323 2313 2312 1313 1312 1212'
  //here we use index 1,2,3
  tangent[0]  = dSdE(1,1,1,1); tangent[1]  = dSdE(1,1,2,2); tangent[2]  = dSdE(1,1,3,3);
  tangent[3]  = dSdE(1,1,2,3); tangent[4]  = dSdE(1,1,1,3); tangent[5]  = dSdE(1,1,1,2);
  tangent[6]  = dSdE(2,2,2,2); tangent[7]  = dSdE(2,2,3,3); tangent[8]  = dSdE(2,2,2,3);
  tangent[9]  = dSdE(2,2,1,3); tangent[10] = dSdE(2,2,1,2); tangent[11] = dSdE(3,3,3,3);
  tangent[12] = dSdE(3,3,2,3); tangent[13] = dSdE(3,3,1,3); tangent[14] = dSdE(3,3,1,2);
  tangent[15] = dSdE(2,3,2,3); tangent[16] = dSdE(2,3,1,3); tangent[17] = dSdE(2,3,1,2);
  tangent[18] = dSdE(1,3,1,3); tangent[19] = dSdE(1,3,1,2); tangent[20] = dSdE(1,2,1,2);
}