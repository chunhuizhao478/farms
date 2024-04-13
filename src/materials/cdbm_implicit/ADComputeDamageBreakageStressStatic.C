//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADComputeDamageBreakageStressStatic.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include "SymmetricRankTwoTensor.h"
#include "SymmetricRankFourTensor.h"

#include <random>

registerMooseObject("farmsApp", ADComputeDamageBreakageStressStatic);
//registerMooseObject("farmsApp", ADSymmetricDamageBreakageStress); //?

//template <typename RankTwoTensor, typename R4, typename ADReal>
InputParameters
ADComputeDamageBreakageStressStatic::validParams()
{
  InputParameters params = ADComputeDamageBreakageStressBase::validParams();
  params.addClassDescription("Compute stress using elasticity for small strains");

  //constant parameters
  params.addRequiredParam<Real>(            "xi_0", "strain invariants ratio: onset of damage evolution");
  params.addRequiredParam<Real>(            "xi_d", "strain invariants ratio: onset of breakage healing");
  params.addRequiredParam<Real>(            "xi_1", "critical point of three phases");
  params.addRequiredParam<Real>(          "xi_min", "strain invariants ratio: minimum allowable value");
  params.addRequiredParam<Real>(          "xi_max", "strain invariants ratio: maximum allowable value");
  params.addRequiredParam<Real>(              "a0", "parameters in granular states");
  params.addRequiredParam<Real>(              "a1", "parameters in granular states");
  params.addRequiredParam<Real>(              "a2", "parameters in granular states");
  params.addRequiredParam<Real>(              "a3", "parameters in granular states");
  params.addRequiredParam<Real>( "gamma_damaged_r", "coefficient of damage solid modulus");
  params.addRequiredParam<Real>(             "C_g", "material parameter: compliance or fluidity of the fine grain granular material");
  params.addRequiredParam<Real>(              "m1", "coefficient of std::power law indexes");
  params.addRequiredParam<Real>(              "m2", "coefficient of std::power law indexes");
  params.addRequiredParam<Real>(               "D", "coefficient of alpha gradient");

  //variable parameters
  params.addRequiredCoupledVar("alpha_in", "damage variable computed from subApp");
  params.addRequiredCoupledVar(    "B_in", "breakage variable computed from subApp");
  params.addRequiredCoupledVar("alpha_grad_x", "damage variable gradient component in x computed from subApp");
  params.addRequiredCoupledVar("alpha_grad_y", "damage variable gradient component in y computed from subApp");
  params.addRequiredCoupledVar("alpha_grad_z", "damage variable gradient component in z computed from subApp");
  params.addRequiredCoupledVar("initial_alpha", "initial distribution of alpha");  

  return params;
}

//template <typename RankTwoTensor, typename R4, typename ADReal>
ADComputeDamageBreakageStressStatic::ADComputeDamageBreakageStressStatic(
    const InputParameters & parameters)
  : ADComputeDamageBreakageStressBase(parameters),
    //_elasticity_tensor_name(_base_name + "elasticity_tensor"),
    //_elasticity_tensor(getADMaterialProperty<RankFourTensor>(_elasticity_tensor_name)),
    _static_initial_stress_tensor(getADMaterialProperty<RankTwoTensor>("static_initial_stress_tensor")),
    _static_initial_strain_tensor(getADMaterialProperty<RankTwoTensor>("static_initial_strain_tensor")),
    _xi_0(getParam<Real>("xi_0")),
    _xi_d(getParam<Real>("xi_d")),
    _xi_1(getParam<Real>("xi_1")),
    _xi_min(getParam<Real>("xi_min")),
    _xi_max(getParam<Real>("xi_max")),
    _a0(getParam<Real>("a0")),
    _a1(getParam<Real>("a1")),
    _a2(getParam<Real>("a2")),
    _a3(getParam<Real>("a3")),
    _gamma_damaged_r(getParam<Real>("gamma_damaged_r")),
    _C_g(getParam<Real>("C_g")),
    _m1(getParam<Real>("m1")),
    _m2(getParam<Real>("m2")),
    _alpha_damagedvar_old(getMaterialPropertyOld<Real>("alpha_damagedvar")),
    _B_old(getMaterialPropertyOld<Real>("B")),
    _xi_old(getMaterialPropertyOld<Real>("xi")),
    _I1_old(getMaterialPropertyOld<Real>("I1")),
    _I2_old(getMaterialPropertyOld<Real>("I2")),
    _lambda_old(getMaterialPropertyOld<Real>("lambda")),
    _shear_modulus_old(getMaterialPropertyOld<Real>("shear_modulus")),
    _gamma_damaged_old(getMaterialPropertyOld<Real>("gamma_damaged")),
    _eps_total_old(getMaterialPropertyOld<RankTwoTensor>("eps_total")),
    _mechanical_strain_old(getMaterialPropertyOld<RankTwoTensor>("mechanical_strain")),
    _eps_p_old(getMaterialPropertyOld<RankTwoTensor>("eps_p")),
    _eps_e_old(getMaterialPropertyOld<RankTwoTensor>("eps_e")),
    _alpha_in(adCoupledValue("alpha_in")),
    _B_in(adCoupledValue("B_in")),
    _alpha_grad_x(adCoupledValue("alpha_grad_x")),
    _alpha_grad_y(adCoupledValue("alpha_grad_y")),
    _alpha_grad_z(adCoupledValue("alpha_grad_z")),
    _density_old(getADMaterialProperty<Real>("density")),
    _D(getParam<Real>("D")),
    _initial_alpha(adCoupledValue("initial_alpha"))      
{
}

//template <typename RankTwoTensor, typename R4, typename ADReal>
void
ADComputeDamageBreakageStressStatic::initialSetup()
{
  if ( hasBlockMaterialProperty<RankTwoTensor>(_base_name + "strain_increment"))
    mooseError("This linear elastic stress calculation only works for small strains; use "
               "ADComputeFiniteStrainElasticStress for simulations using incremental and finite "
               "strains.");
}

//template <typename RankTwoTensor, typename R4, typename ADReal>
void
ADComputeDamageBreakageStressStatic::computeQpStress()
{

  /* 
  compute alpha and B parameters
  */
  ADReal x_coord = _q_point[_qp](0); //along the strike direction
  ADReal y_coord = _q_point[_qp](1); //along the normal direction
  
  ADReal alpha_o; 
  std::random_device rd;
  std::mt19937 gen(rd());
  std::weibull_distribution<double> wb_distribution(2.0,0.05);

  if (y_coord >= 0-1*0.1 and y_coord <= 0+1*0.1){
    alpha_o = 0.7;
  }
  else if (y_coord >= -45 and y_coord <= 45){
    alpha_o = wb_distribution(gen);
  }
  else{
    alpha_o = 0.0;
  }

  //alpha, B are updated
  ADReal alpha_out = alpha_o;
  ADReal B_out = 0; //no breakage evolution

  //grad_alpha
  ADReal alpha_grad_x = _alpha_grad_x[_qp];
  ADReal alpha_grad_y = _alpha_grad_y[_qp];
  ADReal alpha_grad_z = _alpha_grad_z[_qp];
  ADReal D = _D;

  //save alpha and B
  _alpha_damagedvar[_qp] = _alpha_in[_qp];
  _B[_qp] = _B_in[_qp];

  /*
  update modulus
  */

  //lambda, shear_modulus, gamma_damaged are updated
  ADReal lambda_out = _lambda_o;
  ADReal shear_modulus_out = _shear_modulus_o + alpha_out * _xi_0 * _gamma_damaged_r;
  ADReal gamma_damaged_out = alpha_out * _gamma_damaged_r;

  //save
  _lambda[_qp] = lambda_out;
  _shear_modulus[_qp] = shear_modulus_out;
  _gamma_damaged[_qp] = gamma_damaged_out;  

  /*
    compute strain
  */
    
  //take updated parameters
  ADReal lambda = lambda_out;
  ADReal shear_modulus = shear_modulus_out;
  ADReal gamma_damaged = gamma_damaged_out;
  ADReal B = B_out;

  //Define components
  ADReal eps11t = _mechanical_strain[_qp](0,0);
  ADReal eps22t = _mechanical_strain[_qp](1,1);
  ADReal eps33t = _mechanical_strain[_qp](2,2);
  ADReal eps12t = _mechanical_strain[_qp](0,1);
  ADReal eps13t = _mechanical_strain[_qp](0,2);
  ADReal eps23t = _mechanical_strain[_qp](1,2);

  //compute elastic strain
  ADReal eps11e = eps11t;
  ADReal eps22e = eps22t;
  ADReal eps33e = eps33t;
  ADReal eps12e = eps12t;
  ADReal eps13e = eps13t;
  ADReal eps23e = eps23t;

  //represent I1 I2 xi using elastic strain
  ADReal I1 = eps11e + eps22e + eps33e;
  ADReal I2 = eps11e * eps11e + eps22e * eps22e + eps33e * eps33e + 2 * eps12e * eps12e + 2 * eps13e * eps13e + 2 * eps23e * eps23e;
  ADReal xi = I1 / std::sqrt(I2);

  //Represent sigma (solid(s) + granular(b))
  ADReal sigma11_s = ( lambda - gamma_damaged / xi ) * I1 + ( 2 * shear_modulus - gamma_damaged * xi ) * eps11e;
  ADReal sigma11_b = ( 2 * _a2 + _a1 / xi + 3 * _a3 * xi ) * I1 + ( 2 * _a0 + _a1 * xi - _a3 * std::pow(xi,3) ) * eps11e;
  ADReal sigma22_s = ( lambda - gamma_damaged / xi ) * I1 + ( 2 * shear_modulus - gamma_damaged * xi ) * eps22e;
  ADReal sigma22_b = ( 2 * _a2 + _a1 / xi + 3 * _a3 * xi ) * I1 + ( 2 * _a0 + _a1 * xi - _a3 * std::pow(xi,3) ) * eps22e;
  ADReal sigma33_s = ( lambda - gamma_damaged / xi ) * I1 + ( 2 * shear_modulus - gamma_damaged * xi ) * eps33e;
  ADReal sigma33_b = ( 2 * _a2 + _a1 / xi + 3 * _a3 * xi ) * I1 + ( 2 * _a0 + _a1 * xi - _a3 * std::pow(xi,3) ) * eps33e;
  ADReal sigma12_s = ( 2 * shear_modulus - gamma_damaged * xi ) * eps12e;
  ADReal sigma12_b = ( 2 * _a0 + _a1 * xi - _a3 * std::pow(xi,3) ) * eps12e;
  ADReal sigma13_s = ( 2 * shear_modulus - gamma_damaged * xi ) * eps13e;
  ADReal sigma13_b = ( 2 * _a0 + _a1 * xi - _a3 * std::pow(xi,3) ) * eps13e;
  ADReal sigma23_s = ( 2 * shear_modulus - gamma_damaged * xi ) * eps23e;
  ADReal sigma23_b = ( 2 * _a0 + _a1 * xi - _a3 * std::pow(xi,3) ) * eps23e;

  //Represent total stress
  ADReal sigma11_t = (1 - B) * sigma11_s + B * sigma11_b;
  ADReal sigma22_t = (1 - B) * sigma22_s + B * sigma22_b;
  ADReal sigma33_t = (1 - B) * sigma33_s + B * sigma33_b;
  ADReal sigma12_t = (1 - B) * sigma12_s + B * sigma12_b;
  ADReal sigma13_t = (1 - B) * sigma13_s + B * sigma13_b;
  ADReal sigma23_t = (1 - B) * sigma23_s + B * sigma23_b;

  //save data
  //elastic strain
  ADRankTwoTensor eps_e_out;
  ADReal eps11e_out = eps11t;
  ADReal eps22e_out = eps22t;
  ADReal eps33e_out = eps33t;
  ADReal eps12e_out = eps12t;
  ADReal eps13e_out = eps13t;
  ADReal eps23e_out = eps23t;
  eps_e_out(0,0) = eps11e_out;
  eps_e_out(1,1) = eps22e_out;
  eps_e_out(0,1) = eps12e_out;
  eps_e_out(1,0) = eps12e_out;
  eps_e_out(0,2) = eps13e_out;
  eps_e_out(2,0) = eps13e_out;
  eps_e_out(1,2) = eps23e_out;
  eps_e_out(2,1) = eps23e_out;
  eps_e_out(2,2) = eps33e_out;
  _eps_e[_qp] = eps_e_out;

  //total strain
  ADRankTwoTensor eps_total_out;
  eps_total_out(0,0) = eps11t;
  eps_total_out(1,1) = eps22t;
  eps_total_out(0,1) = eps12t;
  eps_total_out(1,0) = eps12t;
  eps_total_out(0,2) = eps13t;
  eps_total_out(2,0) = eps13t;
  eps_total_out(1,2) = eps23t;
  eps_total_out(2,1) = eps23t;
  eps_total_out(2,2) = eps33t;
  _eps_total[_qp] = eps_total_out;

  /*
    compute I1 I2 xi
  */

  //save
  _I1[_qp] = I1;
  _I2[_qp] = I2;
  _xi[_qp] = xi; 

  //feed stress change (relative to initial condition) to system
  ADRankTwoTensor stress_out;
  stress_out(0,0) = sigma11_t;
  stress_out(1,1) = sigma22_t;
  stress_out(2,2) = sigma33_t;
  stress_out(0,1) = sigma12_t;
  stress_out(1,0) = sigma12_t;
  stress_out(0,2) = sigma13_t;
  stress_out(2,0) = sigma13_t;
  stress_out(1,2) = sigma23_t;
  stress_out(2,1) = sigma23_t; 

  // stress = C * e
  _stress[_qp] = stress_out;

  // Assign value for elastic strain, which is equal to the mechanical strain
  _elastic_strain[_qp] = _mechanical_strain[_qp];
}

//template class ADComputeDamageBreakageStressTempl<RankTwoTensor, RankFourTensor, ADReal>;
//template class ADComputeDamageBreakageStressTempl<SymmetricRankTwoTensor, SymmetricRankFourTensor>;