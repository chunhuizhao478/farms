//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ComputeDamageBreakageStress3DDynamicCDBM.h"
#include "NestedSolve.h"
#include "FEProblem.h"

registerMooseObject("farmsApp", ComputeDamageBreakageStress3DDynamicCDBM);

InputParameters
ComputeDamageBreakageStress3DDynamicCDBM::validParams()
{ 
  //Note: lambda_o, shear_modulus_o is defined in "ComputeGeneralDamageBreakageStressBase"
  //to initialize _lambda, _shear_modulus material properties
  InputParameters params = ComputeDamageBreakageStressBase3D::validParams();
  params.addClassDescription("Compute stress using elasticity for small strains");
  
  //constant parameters
  params.addRequiredParam<Real>(        "lambda_o", "initial lambda constant value");
  params.addRequiredParam<Real>( "shear_modulus_o", "initial shear modulus value");
  params.addRequiredParam<Real>(            "xi_0", "strain invariants ratio: onset of damage evolution");
  params.addRequiredParam<Real>(            "xi_d", "strain invariants ratio: onset of breakage healing");
  params.addRequiredParam<Real>(          "xi_min", "strain invariants ratio: minimum allowable value");
  params.addRequiredParam<Real>(          "xi_max", "strain invariants ratio: maximum allowable value");
  params.addRequiredParam<Real>(             "chi", "ratio of solid energy and granular energy");
  params.addRequiredParam<Real>(             "C_g", "material parameter: compliance or fluidity of the fine grain granular material");
  params.addRequiredParam<Real>(              "m1", "coefficient of std::power law indexes");
  params.addRequiredParam<Real>(              "m2", "coefficient of std::power law indexes");
  params.addRequiredParam<Real>(     "Cd_constant", "coefficient gives positive damage evolution");
  params.addRequiredParam<Real>(             "C_1", "coefficient of healing for damage evolution");
  params.addRequiredParam<Real>(             "C_2", "coefficient of healing for damage evolution");
  params.addRequiredParam<Real>(      "beta_width", "coefficient gives width of transitional region");
  params.addRequiredParam<Real>( "CdCb_multiplier", "multiplier between Cd and Cb");
  params.addRequiredParam<Real>(    "CBH_constant", "constant CBH value");
  params.addRequiredParam<Real>(    "D", "D value");
  
  //variable parameters
  params.addRequiredCoupledVar("alpha_grad_x", "damage variable gradient component in x computed from subApp");
  params.addRequiredCoupledVar("alpha_grad_y", "damage variable gradient component in y computed from subApp");
  params.addRequiredCoupledVar("alpha_grad_z", "damage variable gradient component in z computed from subApp");

  return params;
}

ComputeDamageBreakageStress3DDynamicCDBM::ComputeDamageBreakageStress3DDynamicCDBM(const InputParameters & parameters)
  : ComputeDamageBreakageStressBase3D(parameters),
    _xi_0(getParam<Real>("xi_0")),
    _xi_d(getParam<Real>("xi_d")),
    _xi_min(getParam<Real>("xi_min")),
    _xi_max(getParam<Real>("xi_max")),
    _chi(getParam<Real>("chi")),
    _C_g(getParam<Real>("C_g")),
    _m1(getParam<Real>("m1")),
    _m2(getParam<Real>("m2")),
    _alpha_damagedvar_old(getMaterialPropertyOldByName<Real>("alpha_damagedvar")),
    _B_old(getMaterialPropertyOldByName<Real>("B")),
    _xi_old(getMaterialPropertyOldByName<Real>("xi")),
    _I1_old(getMaterialPropertyOldByName<Real>("I1")),
    _I2_old(getMaterialPropertyOldByName<Real>("I2")),
    _lambda_old(getMaterialPropertyOldByName<Real>("lambda")),
    _shear_modulus_old(getMaterialPropertyOldByName<Real>("shear_modulus")),
    _gamma_damaged_old(getMaterialPropertyOldByName<Real>("gamma_damaged")),
    _eps_total_old(getMaterialPropertyOldByName<RankTwoTensor>("eps_total")),
    _mechanical_strain_old(getMaterialPropertyOldByName<RankTwoTensor>("mechanical_strain")),
    _eps_p_old(getMaterialPropertyOldByName<RankTwoTensor>("eps_p")),
    _eps_e_old(getMaterialPropertyOldByName<RankTwoTensor>("eps_e")),
    _sigma_d_old(getMaterialPropertyOldByName<RankTwoTensor>("sigma_d")),
    _alpha_grad_x(coupledValue("alpha_grad_x")),
    _alpha_grad_y(coupledValue("alpha_grad_y")),
    _alpha_grad_z(coupledValue("alpha_grad_z")),
    _D(getParam<Real>("D")),
    _initial_damage(getMaterialPropertyByName<Real>("initial_damage")),
    _initial_breakage(getMaterialPropertyByName<Real>("initial_breakage")),
    _initial_shear_stress(getMaterialPropertyByName<Real>("initial_shear_stress")),
    _damage_perturbation(getMaterialPropertyByName<Real>("damage_perturbation")),
    _shear_stress_perturbation(getMaterialPropertyByName<Real>("shear_stress_perturbation")),
    _Cd_constant(getParam<Real>("Cd_constant")),
    _C1(getParam<Real>("C_1")),
    _C2(getParam<Real>("C_2")),
    _beta_width(getParam<Real>("beta_width")),
    _CdCb_multiplier(getParam<Real>("CdCb_multiplier")),
    _CBH_constant(getParam<Real>("CBH_constant")),
    _dim(_mesh.dimension())
{
}

void
ComputeDamageBreakageStress3DDynamicCDBM::initialSetup()
{
  // _base_name + "unstabilized_deformation_gradient" is only declared if we're
  // using the Lagrangian kernels.  It's okay to invoke this small strain
  // material if you are using that kernel system and the
  // ComputeLagrangianWrappedStress wrapper
  if (hasBlockMaterialProperty<RankTwoTensor>(_base_name + "strain_increment") &&
      !hasBlockMaterialProperty<RankTwoTensor>(_base_name + "unstabilized_deformation_gradient"))
    mooseError("This linear elastic stress calculation only works for small strains; use "
               "ComputeFiniteStrainElasticStress for simulations using incremental and finite "
               "strains.");
               
}

void
ComputeDamageBreakageStress3DDynamicCDBM::initQpStatefulProperties()
{
  _elastic_strain[_qp].zero();
  _stress[_qp].zero();
  // _alpha_damagedvar[_qp] = _initial_damage[_qp];

}

void
ComputeDamageBreakageStress3DDynamicCDBM::computeQpStress()
{ 
  
  /*
  compute gammar, breakage coefficients
  */
  Real gamma_damaged_r = computegammar();
  std::vector<Real> avec = computecoefficients(gamma_damaged_r);
  Real a0 = avec[0];
  Real a1 = avec[1];
  Real a2 = avec[2];
  Real a3 = avec[3];

  /* 
  compute alpha and B parameters
  */

  /* compute alpha */
  //compute forcing term
  Real alpha_forcingterm;
  if ( _xi_old[_qp] >= _xi_0 && _xi_old[_qp] <= _xi_max ){
    alpha_forcingterm = (1 - _B_old[_qp]) * ( _Cd_constant * _I2_old[_qp] * ( _xi_old[_qp] - _xi_0 ) );
  }
  else if ( _xi_old[_qp] < _xi_0 && _xi_old[_qp] >= _xi_min ){
    alpha_forcingterm = (1 - _B_old[_qp]) * ( _C1 * std::exp(_alpha_damagedvar_old[_qp]/_C2) * _I2_old[_qp] * ( _xi_old[_qp] - _xi_0 ) );
  }
  else{
    mooseError("xi_old is OUT-OF-RANGE!.");   
  }

  //update alpha at current time
  Real alpha_out = _alpha_damagedvar_old[_qp] + _dt * alpha_forcingterm;

  //check alpha within range
  if ( alpha_out < 0 ){ alpha_out = 0.0; }
  else if ( alpha_out > 1 ){ alpha_out = 1.0; }
  else{}       

  //check below initial damage (fix initial damage)
  if ( alpha_out < _initial_damage[_qp] + _damage_perturbation[_qp] ){ alpha_out = _initial_damage[_qp] + _damage_perturbation[_qp]; }
  else{}

  _alpha_damagedvar[_qp] = alpha_out;

  /* compute B */
  Real C_B = _CdCb_multiplier * _Cd_constant;

  //compute xi_1
  Real _xi_1 = _xi_0 + sqrt( pow(_xi_0 , 2) + 2 * _shear_modulus_o / _lambda_o );

  //alphacr function
  Real alphacr;
  if ( _xi_old[_qp] < _xi_0 ){ alphacr = 1.0;} 
  else if ( _xi_old[_qp] > _xi_0 && _xi_old[_qp] <= _xi_1 ){ alphacr = alphacr_root1(_xi_old[_qp],gamma_damaged_r);}
  else if ( _xi_old[_qp] > _xi_1 && _xi_old[_qp] <= _xi_max ){ alphacr = alphacr_root2(_xi_old[_qp],gamma_damaged_r); }
  else{std::cout<<"xi: "<<_xi_old[_qp]<<std::endl;mooseError("xi exceeds the maximum allowable range!");}

  //compute forcing func
  Real Prob = 1.0 / ( std::exp( (alphacr - _alpha_damagedvar_old[_qp]) / _beta_width ) + 1.0 );
  Real B_forcingterm;
  if ( _xi_old[_qp] >= _xi_d && _xi_old[_qp] <= _xi_max ){
    B_forcingterm = 1.0 * C_B * Prob * (1-_B_old[_qp]) * _I2_old[_qp] * (_xi_old[_qp] - _xi_d); //could heal if xi < xi_0
  }
  else if ( _xi_old[_qp] < _xi_d && _xi_old[_qp] >= _xi_min ){
    B_forcingterm = 1.0 * _CBH_constant * _I2_old[_qp] * ( _xi_old[_qp] - _xi_d ); //close healing
  }
  else{
    mooseError("xi_old is OUT-OF-RANGE!.");
  }

  Real B_out = _B_old[_qp] + _dt * B_forcingterm;

  //check breakage within range
  if ( B_out < 0 ){ B_out = 0.0; }
  else if ( B_out > 1 ){ B_out = 1.0; }
  else{}   

  //check below initial damage (fix initial damage)
  if ( B_out < _initial_breakage[_qp] ){ B_out = _initial_breakage[_qp]; }
  else{}

  //save alpha and B
  _B[_qp] = B_out;

  //lambda, shear_modulus, gamma_damaged are updated
  Real lambda_out = _lambda_o;
  Real shear_modulus_out = _shear_modulus_o + alpha_out * _xi_0 * gamma_damaged_r;
  Real gamma_damaged_out = alpha_out * gamma_damaged_r;

  //save
  _lambda[_qp] = lambda_out;
  _shear_modulus[_qp] = shear_modulus_out;
  _gamma_damaged[_qp] = gamma_damaged_out;

  /* compute strain */
  RankTwoTensor eps_p = _eps_p_old[_qp] + _dt * _C_g * std::pow(_B_old[_qp],_m1) * _sigma_d_old[_qp];
  RankTwoTensor eps_e = _mechanical_strain[_qp] - eps_p;

  const Real epsilon = 1e-12;
  Real I1 = epsilon + eps_e(0,0) + eps_e(1,1) + eps_e(2,2);
  Real I2 = epsilon + eps_e(0,0) * eps_e(0,0) + eps_e(1,1) * eps_e(1,1) + eps_e(2,2) * eps_e(2,2) + 2 * eps_e(0,1) * eps_e(0,1) + 2 * eps_e(0,2) * eps_e(0,2) + 2 * eps_e(1,2) * eps_e(1,2);
  Real xi = I1/std::sqrt(I2);

  //Represent sigma (solid(s) + granular(b))
  RankTwoTensor sigma_s;
  RankTwoTensor sigma_b;
  RankTwoTensor sigma_total;
  RankTwoTensor sigma_d;
  const auto I = RankTwoTensor::Identity();

  /* Compute stress */
  sigma_s = (lambda_out - gamma_damaged_out / xi) * I1 * RankTwoTensor::Identity() + (2 * shear_modulus_out - gamma_damaged_out * xi) * eps_e;
  sigma_b = (2 * a2 + a1 / xi + 3 * a3 * xi) * I1 * RankTwoTensor::Identity() + (2 * a0 + a1 * xi - a3 * std::pow(xi, 3)) * eps_e;
  sigma_total = (1 - B_out) * sigma_s + B_out * sigma_b;

  //Add shear perturbation
  if (_shear_stress_perturbation[_qp] != 0){
    sigma_total(0,2) = _initial_shear_stress[_qp] + _shear_stress_perturbation[_qp];
    sigma_total(2,0) = _initial_shear_stress[_qp] + _shear_stress_perturbation[_qp];
  }

  sigma_d = sigma_total - 1/3 * (sigma_total(0,0) + sigma_total(1,1) + sigma_total(2,2)) * I;

  _eps_total[_qp] = eps_p + eps_e;
  _eps_p[_qp] = eps_p;
  _eps_e[_qp] = eps_e;
  _I1[_qp] = I1;
  _I2[_qp] = I2;
  _xi[_qp] = xi;
  _sigma_d[_qp] = sigma_d;

  // Rotate the stress state to the current configuration
  _stress[_qp] = sigma_total;

  // Assign value for elastic strain, which is equal to the mechanical strain
  _elastic_strain[_qp] = eps_e;

  // Compute jacobian //_Jacobian_mult[_qp]
  // computeQpTangentModulus(I1, 
  //                         I2, 
  //                         xi, 
  //                         B_out,
  //                         shear_modulus_out, 
  //                         gamma_damaged_out, 
  //                         a0, a1, a2, a3, eps_e);

  //Compute equivalent strain rate
  RankTwoTensor epsilon_rate = (eps_p - _eps_p_old[_qp])/_dt;
  Real epsilon_eq = sqrt(2/3*(epsilon_rate(0,0)*epsilon_rate(0,0)+epsilon_rate(1,1)*epsilon_rate(1,1)+epsilon_rate(2,2)*epsilon_rate(2,2)+2*epsilon_rate(0,1)*epsilon_rate(0,1)+2*epsilon_rate(0,2)*epsilon_rate(0,2)+2*epsilon_rate(1,2)*epsilon_rate(1,2)));
  _epsilon_eq[_qp] = epsilon_eq;

}

Real 
ComputeDamageBreakageStress3DDynamicCDBM::computegammar()
{
  // Calculate each part of the expression
  Real term1 = -_xi_0 * (-_lambda_o * pow(_xi_0, 2) + 6 * _lambda_o + 2 * _shear_modulus_o);
  Real term2_sqrt = sqrt((_lambda_o * pow(_xi_0, 2) + 2 * _shear_modulus_o) * 
                            (_lambda_o * pow(_xi_0, 4) - 12 * _lambda_o * pow(_xi_0, 2) + 36 * _lambda_o 
                            - 6 * _shear_modulus_o * pow(_xi_0, 2) + 24 * _shear_modulus_o));
  Real denominator = 2 * (pow(_xi_0, 2) - 3);
  
  // Calculate gamma_r
  Real gamma_r = (term1 - term2_sqrt) / denominator;
  
  //save
  return gamma_r;
}

std::vector<Real>
ComputeDamageBreakageStress3DDynamicCDBM::computecoefficients(Real gamma_damaged_r)
{

  //compute xi_1
  Real _xi_1 = _xi_0 + sqrt( pow(_xi_0 , 2) + 2 * _shear_modulus_o / _lambda_o );

  //compute alpha_cr | xi = 0
  Real alpha_cr_xi0 = alphacr_root1(0, gamma_damaged_r);

  //compute mu_cr
  Real mu_cr = _shear_modulus_o + alpha_cr_xi0 * _xi_0 * gamma_damaged_r;

  //a0
  Real a0 = _chi * mu_cr;

  //a1
  Real numerator_a1 = -2 * _chi * mu_cr * pow(_xi_1, 3) + 6 * _chi * mu_cr * _xi_1 * pow(_xi_d, 2) - 4 * _chi * mu_cr * pow(_xi_d, 3)
                      - 2 * gamma_damaged_r * pow(_xi_1, 3) * _xi_d + 2 * gamma_damaged_r * pow(_xi_1, 3) * _xi_0
                      + _lambda_o * pow(_xi_1, 3) * pow(_xi_d, 2) + 2 * _shear_modulus_o * pow(_xi_1, 3);
  Real denominator_a1 = 2 * pow(_xi_1, 3) * _xi_d - 4 * pow(_xi_1, 2) * pow(_xi_d, 2) + 2 * _xi_1 * pow(_xi_d, 3);
  Real a1 = numerator_a1 / denominator_a1;

  //a2
  Real numerator_a2 = 2 * _chi * mu_cr * pow(_xi_1, 3) - 3 * _chi * mu_cr * pow(_xi_1, 2) * _xi_d + _chi * mu_cr * pow(_xi_d, 3)
                       + 2 * gamma_damaged_r * pow(_xi_1, 3) * _xi_d - 2 * gamma_damaged_r * pow(_xi_1, 3) * _xi_0
                       - _lambda_o * pow(_xi_1, 3) * pow(_xi_d, 2) - 2 * _shear_modulus_o * pow(_xi_1, 3);
  Real denominator_a2 = pow(_xi_1, 4) * _xi_d - 2 * pow(_xi_1, 3) * pow(_xi_d, 2) + pow(_xi_1, 2) * pow(_xi_d, 3); 
  Real a2 = numerator_a2 / denominator_a2; 

  //a3
  Real numerator_a3 = -2 * _chi * mu_cr * pow(_xi_1, 2) + 4 * _chi * mu_cr * _xi_1 * _xi_d - 2 * _chi * mu_cr * pow(_xi_d, 2)
                       - 2 * gamma_damaged_r * pow(_xi_1, 2) * _xi_d + 2 * gamma_damaged_r * pow(_xi_1, 2) * _xi_0
                       + _lambda_o * pow(_xi_1, 2) * pow(_xi_d, 2) + 2 * _shear_modulus_o * pow(_xi_1, 2);
  Real denominator_a3 = 2 * pow(_xi_1, 4) * _xi_d - 4 * pow(_xi_1, 3) * pow(_xi_d, 2) + 2 * pow(_xi_1, 2) * pow(_xi_d, 3);
  Real a3 = numerator_a3 / denominator_a3; 

  //save
  std::vector<Real> a_vec {a0,a1,a2,a3};

  return a_vec;

}

// Function for alpha_func_root1
Real 
ComputeDamageBreakageStress3DDynamicCDBM::alphacr_root1(Real xi, Real gamma_damaged_r) {
    Real term1 = _lambda_o * pow(xi, 3) - 6 * _lambda_o * _xi_0 + 6 * _shear_modulus_o * xi - 8 * _shear_modulus_o * _xi_0;
    Real term2 = std::sqrt(_lambda_o * _lambda_o * pow(xi, 6) 
                             - 12 * _lambda_o * _lambda_o * pow(xi, 3) * _xi_0 
                             + 36 * _lambda_o * _lambda_o * _xi_0 * _xi_0 
                             + 12 * _lambda_o * _shear_modulus_o * pow(xi, 4) 
                             - 16 * _lambda_o * _shear_modulus_o * pow(xi, 3) * _xi_0 
                             - 72 * _lambda_o * _shear_modulus_o * pow(xi, 2) 
                             + 72 * _lambda_o * _shear_modulus_o * xi * _xi_0 
                             + 72 * _lambda_o * _shear_modulus_o 
                             - 12 * _shear_modulus_o * _shear_modulus_o * pow(xi, 2) 
                             + 48 * _shear_modulus_o * _shear_modulus_o);
    Real denominator = 2 * gamma_damaged_r * (3 * pow(xi, 2) - 6 * xi * _xi_0 + 4 * _xi_0 * _xi_0 - 3);
    return (term1 - term2) / denominator;
}

// Function for alpha_func_root2
Real 
ComputeDamageBreakageStress3DDynamicCDBM::alphacr_root2(Real xi, Real gamma_damaged_r) {
    return 2 * _shear_modulus_o / (gamma_damaged_r * (xi - 2 * _xi_0));
}

// void
// ComputeDamageBreakageStress3DDynamicCDBM::computeQpTangentModulus(Real I1, 
//                                                        Real I2, 
//                                                        Real xi, 
//                                                        Real B,
//                                                        Real shear_modulus_out,
//                                                        Real gamma_damaged_out,
//                                                        Real a0,
//                                                        Real a1,
//                                                        Real a2,
//                                                        Real a3,
//                                                        RankTwoTensor Ee)
// {

//   // //define functions for derivatives dstress_dstrain
//   RankFourTensor tangent;

//   //delta function
//   auto delta = [](int i, int j) -> Real {
//     return (i == j) ? 1.0 : 0.0;
//   };

//   //dI1_dE_{kl}
//   auto dI1dE = [&](int k, int l) -> Real {
//     return delta(k,l);
//   };

//   //dI2_dE_{kl}
//   auto dI2dE = [&](int k, int l) -> Real {
//     return 2 * Ee(k,l);
//   };

//   //dxi_dE_{kl}
//   auto dxidE = [&](int k, int l) -> Real {
    
//     // Epsilon to avoid division by zero
//     const Real epsilon = 1e-12;
//     // Adjust I2 if necessary
//     Real adjusted_I2 = I2;
//     if (I2 <= epsilon) {
//       //mooseWarning("I2 is zero or too small (I2 = ", I2, "), adjusting to epsilon.");
//       adjusted_I2 = epsilon;
//     }

//     Real dxidE = 0.5 * pow(adjusted_I2,-1.5) * dI2dE(k,l) * I1;
//     //mooseInfo("I1 = ", I1, ", I2 = ", I2);
//     if (std::isnan(dxidE)){mooseError("dxidE");}
//     return delta(k,l) * pow(adjusted_I2,-0.5) - 0.5 * pow(adjusted_I2,-1.5) * dI2dE(k,l) * I1;
//   };

//   //dE_{ij}_dE_{kl}
//   auto dEdE = [&](int i, int j, int k, int l) -> Real {
//     return delta(i,k) * delta(j,l);
//     //return 0.5 * ( delta(i,k) * delta(j,l) + delta(i,l) * delta(j,k) ); //its symmetric form
//   };

//   //dxi^{-1}_dE_{kl}
//   auto dxim1dE = [&](int k, int l) -> Real {
//     return -1.0 * pow(xi,-2.0) * dxidE(k,l);
//   };

//   //dxi^3_dE_{kl}
//   auto dxi3dE = [&](int k, int l) -> Real {
//     return 3 * pow(xi,2) * dxidE(k,l);
//   };

//   //dSe_{ij}_dE_{kl}
//   auto dSedE = [&](int i, int j, int k, int l) -> Real {
//     Real dSedE_components = (- gamma_damaged_out * dxim1dE(k,l) ) * I1 * delta(i,j);
//     dSedE_components += ( _lambda_o - gamma_damaged_out / xi ) * dI1dE(k,l) * delta(i,j);
//     dSedE_components += (- gamma_damaged_out * dxidE(k,l) ) * Ee(i,j);
//     dSedE_components += ( 2 * shear_modulus_out - gamma_damaged_out * xi ) * dEdE(i,j,k,l);
//     if (std::isnan(dSedE_components)){mooseError("dSedE_components");}
//     return dSedE_components;
//   };

//   //dSb_{ij}_dE_{kl}
//   auto dSbdE = [&](int i, int j, int k, int l) -> Real {
//     Real dSbdE_components = ( a1 * dxim1dE(k,l) + 3 * a3 * dxidE(k,l) ) * I1 * delta(i,j);
//     dSbdE_components += ( 2 * a2 + a1 / xi + 3 * a3 * xi ) * dI1dE(k,l) * delta(i,j);
//     dSbdE_components += ( a1 * dxidE(k,l) - a3 * dxi3dE(k,l) ) * Ee(i,j);
//     dSbdE_components += ( 2 * a0 + a1 * xi - a3 * pow(xi,3) ) * dEdE(i,j,k,l);
//     return dSbdE_components;
//   };

//   //dS_{ij}_dE_{kl}
//   auto dSdE = [&](int i, int j, int k, int l) -> Real {
//     return (1 - B) * dSedE(i,j,k,l) + B * dSbdE(i,j,k,l);
//   };

//   // Compute tangent modulus C
//   for (unsigned int i = 0; i < _dim; i++){
//     for (unsigned int j = 0; j < _dim; j++){
//       for (unsigned int k = 0; k < _dim; k++){
//         for (unsigned int l = 0; l < _dim; l++){
//           if (std::isnan(dSdE(i,j,k,l))){mooseError("encounter nan error: dSdE(i,j,k,l)");}
//           tangent(i,j,k,l) += dSdE(i,j,k,l);
//         }
//       }
//     }
//   }

//   //update jacobian
//   _Jacobian_mult[_qp] = tangent;

// }