//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADComputeDamageBreakageStressv2.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include "SymmetricRankTwoTensor.h"
#include "SymmetricRankFourTensor.h"

registerMooseObject("farmsApp", ADComputeDamageBreakageStressv2);
//registerMooseObject("farmsApp", ADSymmetricDamageBreakageStress); //?

//template <typename RankTwoTensor, typename R4, typename ADReal>
InputParameters
ADComputeDamageBreakageStressv2::validParams()
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
ADComputeDamageBreakageStressv2::ADComputeDamageBreakageStressv2(
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
ADComputeDamageBreakageStressv2::initialSetup()
{
  if ( hasBlockMaterialProperty<RankTwoTensor>(_base_name + "strain_increment"))
    mooseError("This linear elastic stress calculation only works for small strains; use "
               "ADComputeFiniteStrainElasticStress for simulations using incremental and finite "
               "strains.");
}

//template <typename RankTwoTensor, typename R4, typename ADReal>
void
ADComputeDamageBreakageStressv2::computeQpStress()
{

  /* 
  compute alpha and B parameters
  */

  if (_t == 0.0)
  {
    setupInitial(); 
    //initialize lambda, shear_modulus, gamma_damaged
    //alpha_damagedvar, B, eps_p, eps_e, eps_total, sts_total, I1, I2, xi,
    //shear_wave_speed, pressure_wave_speed
  }
  else{

    //alpha, B are updated
    ADReal alpha_out = _alpha_in[_qp];
    ADReal B_out = _B_in[_qp];

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

    //eps_p, eps_e, eps_total are updated
    //Retrieve parameter values
    //strain: total, viscoelastic
    ADRankTwoTensor eps_total_change = _mechanical_strain[_qp];
    RankTwoTensor eps_total_change_old = _mechanical_strain_old[_qp];
    ADRankTwoTensor eps_total_change_inc = eps_total_change - eps_total_change_old;

    //get total strain (add increment)
    ADRankTwoTensor eps_total = _eps_total_old[_qp] + eps_total_change_inc;

    //take the viscoelastic strain from the previous time step
    RankTwoTensor eps_viscoelastic = _eps_p_old[_qp];
      
    //take updated parameters
    ADReal lambda = lambda_out;
    ADReal shear_modulus = shear_modulus_out;
    ADReal gamma_damaged = gamma_damaged_out;
    ADReal B = B_out;

    //Define components
    ADReal eps11p_pre = eps_viscoelastic(0,0);
    ADReal eps22p_pre = eps_viscoelastic(1,1);
    ADReal eps33p_pre = eps_viscoelastic(2,2);
    ADReal eps12p_pre = eps_viscoelastic(0,1);
    ADReal eps13p_pre = eps_viscoelastic(0,2);
    ADReal eps23p_pre = eps_viscoelastic(1,2);
    ADReal eps11t = eps_total(0,0);
    ADReal eps22t = eps_total(1,1);
    ADReal eps33t = eps_total(2,2);
    ADReal eps12t = eps_total(0,1);
    ADReal eps13t = eps_total(0,2);
    ADReal eps23t = eps_total(1,2);

    //compute elastic strain
    ADReal eps11e = eps11t - eps11p_pre;
    ADReal eps22e = eps22t - eps22p_pre;
    ADReal eps33e = eps33t - eps33p_pre;
    ADReal eps12e = eps12t - eps12p_pre;
    ADReal eps13e = eps13t - eps13p_pre;
    ADReal eps23e = eps23t - eps23p_pre;

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

    //Represent deviatoric stress
    ADReal sigma_d11 = sigma11_t - 1/3 * (sigma11_t + sigma22_t + sigma33_t);
    ADReal sigma_d22 = sigma22_t - 1/3 * (sigma11_t + sigma22_t + sigma33_t);
    ADReal sigma_d33 = sigma33_t - 1/3 * (sigma11_t + sigma22_t + sigma33_t);
    ADReal sigma_d12 = sigma12_t;
    ADReal sigma_d13 = sigma13_t;
    ADReal sigma_d23 = sigma23_t;

    //Setup equations
    ///residual
    ADReal eps11p_inc = _dt * _C_g * std::pow(B,_m1) * std::pow(sigma_d11,_m2);
    ADReal eps22p_inc = _dt * _C_g * std::pow(B,_m1) * std::pow(sigma_d22,_m2);
    ADReal eps33p_inc = _dt * _C_g * std::pow(B,_m1) * std::pow(sigma_d33,_m2);
    ADReal eps12p_inc = _dt * _C_g * std::pow(B,_m1) * std::pow(sigma_d12,_m2);
    ADReal eps13p_inc = _dt * _C_g * std::pow(B,_m1) * std::pow(sigma_d13,_m2);
    ADReal eps23p_inc = _dt * _C_g * std::pow(B,_m1) * std::pow(sigma_d23,_m2);

    //save data
    //viscoelastic strain
    ADRankTwoTensor eps_p_out;
    eps_p_out(0,0) = eps11p_pre + eps11p_inc;
    eps_p_out(1,1) = eps22p_pre + eps22p_inc;
    eps_p_out(0,1) = eps12p_pre + eps12p_inc;
    eps_p_out(1,0) = eps12p_pre + eps12p_inc;
    eps_p_out(0,2) = eps13p_pre + eps13p_inc;
    eps_p_out(2,0) = eps13p_pre + eps13p_inc;
    eps_p_out(1,2) = eps23p_pre + eps23p_inc;
    eps_p_out(2,1) = eps23p_pre + eps23p_inc;
    eps_p_out(2,2) = eps33p_pre + eps33p_inc;
    _eps_p[_qp] = eps_p_out;

    //elastic strain
    ADRankTwoTensor eps_e_out;
    ADReal eps11e_out = eps11t - ( eps11p_pre + eps11p_inc );
    ADReal eps22e_out = eps22t - ( eps22p_pre + eps22p_inc );
    ADReal eps33e_out = eps33t - ( eps33p_pre + eps33p_inc );
    ADReal eps12e_out = eps12t - ( eps12p_pre + eps12p_inc );
    ADReal eps13e_out = eps13t - ( eps13p_pre + eps13p_inc );
    ADReal eps23e_out = eps23t - ( eps23p_pre + eps23p_inc );
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

    //I1, I2, xi are updated
    
    ADReal I1_out = eps11e_out + eps22e_out + eps33e_out;
    ADReal I2_out = eps11e_out * eps11e_out + eps22e_out * eps22e_out + eps33e_out * eps33e_out + 2 * eps12e_out * eps12e_out + 2 * eps13e_out * eps13e_out + 2 * eps23e_out * eps23e_out;
    ADReal xi_out = I1_out / std::sqrt(I2_out);

    //save
    _I1[_qp] = I1_out;
    _I2[_qp] = I2_out;
    _xi[_qp] = xi_out; 

    //compute principal strain
    //_principal_strain[_qp] = 0.5 * ( eps11e_out + eps22e_out ) + std::sqrt( std::pow( 0.5 * ( eps11e_out - eps22e_out ) , 2) + std::pow(eps12e_out , 2) );
    _principal_strain[_qp] = I1_out;

    //compute stress
    //sts_total, stress are updated
    //feed total stress
    ADRankTwoTensor stress_total_out;
    stress_total_out(0,0) = ( 1 - B ) * ( ( lambda - gamma_damaged / xi_out ) * I1_out + ( 2 * shear_modulus - gamma_damaged * xi_out ) * eps11e_out ) + B * ( ( 2 * _a2 + _a1 / xi_out + 3 * _a3 * xi_out ) * I1_out + ( 2 * _a0 + _a1 * xi_out - _a3 * std::pow(xi_out,3) ) * eps11e_out );
    stress_total_out(1,1) = ( 1 - B ) * ( ( lambda - gamma_damaged / xi_out ) * I1_out + ( 2 * shear_modulus - gamma_damaged * xi_out ) * eps22e_out ) + B * ( ( 2 * _a2 + _a1 / xi_out + 3 * _a3 * xi_out ) * I1_out + ( 2 * _a0 + _a1 * xi_out - _a3 * std::pow(xi_out,3) ) * eps22e_out );
    stress_total_out(2,2) = ( 1 - B ) * ( ( lambda - gamma_damaged / xi_out ) * I1_out + ( 2 * shear_modulus - gamma_damaged * xi_out ) * eps33e_out ) + B * ( ( 2 * _a2 + _a1 / xi_out + 3 * _a3 * xi_out ) * I1_out + ( 2 * _a0 + _a1 * xi_out - _a3 * std::pow(xi_out,3) ) * eps33e_out );
    stress_total_out(0,1) = ( 1 - B ) * ( ( 2 * shear_modulus - gamma_damaged * xi_out ) * eps12e_out ) + B * ( ( 2 * _a0 + _a1 * xi_out - _a3 * std::pow(xi_out,3) ) * eps12e_out );
    stress_total_out(1,0) = ( 1 - B ) * ( ( 2 * shear_modulus - gamma_damaged * xi_out ) * eps12e_out ) + B * ( ( 2 * _a0 + _a1 * xi_out - _a3 * std::pow(xi_out,3) ) * eps12e_out );
    stress_total_out(0,2) = ( 1 - B ) * ( ( 2 * shear_modulus - gamma_damaged * xi_out ) * eps13e_out ) + B * ( ( 2 * _a0 + _a1 * xi_out - _a3 * std::pow(xi_out,3) ) * eps13e_out );
    stress_total_out(2,0) = ( 1 - B ) * ( ( 2 * shear_modulus - gamma_damaged * xi_out ) * eps13e_out ) + B * ( ( 2 * _a0 + _a1 * xi_out - _a3 * std::pow(xi_out,3) ) * eps13e_out );
    stress_total_out(1,2) = ( 1 - B ) * ( ( 2 * shear_modulus - gamma_damaged * xi_out ) * eps23e_out ) + B * ( ( 2 * _a0 + _a1 * xi_out - _a3 * std::pow(xi_out,3) ) * eps23e_out );
    stress_total_out(2,1) = ( 1 - B ) * ( ( 2 * shear_modulus - gamma_damaged * xi_out ) * eps23e_out ) + B * ( ( 2 * _a0 + _a1 * xi_out - _a3 * std::pow(xi_out,3) ) * eps23e_out );

    //_sts_total[_qp] = stress_total_out;

    //feed stress change (relative to initial condition) to system
    ADRankTwoTensor stress_out;
    ADRankTwoTensor stress_initial = _static_initial_stress_tensor[_qp];
    stress_out(0,0) = stress_total_out(0,0) - stress_initial(0,0);
    stress_out(1,1) = stress_total_out(1,1) - stress_initial(1,1);
    stress_out(2,2) = stress_total_out(2,2) - stress_initial(2,2);
    stress_out(0,1) = stress_total_out(0,1) - stress_initial(0,1);
    stress_out(1,0) = stress_total_out(1,0) - stress_initial(1,0);
    stress_out(0,2) = stress_total_out(0,2) - stress_initial(0,2);
    stress_out(2,0) = stress_total_out(2,0) - stress_initial(2,0);
    stress_out(1,2) = stress_total_out(1,2) - stress_initial(1,2);
    stress_out(2,1) = stress_total_out(2,1) - stress_initial(2,1);

    // stress = C * e
    _stress[_qp] = stress_out;

    // Assign value for elastic strain, which is equal to the mechanical strain
    _elastic_strain[_qp] = _mechanical_strain[_qp];
  }
}

void
ADComputeDamageBreakageStressv2::setupInitial()
{

  ADReal x_coord = _q_point[_qp](0); //along the strike direction
  ADReal y_coord = _q_point[_qp](1); //along the normal direction

  ADReal alpha_o = 0;

  if (y_coord >= 0-2*0.01 and y_coord <= 0+2*0.01){
    if (x_coord >= -0.5-2*0.01 and x_coord <= -0.5+2*0.01){
        alpha_o = 0.8;
    }
    else if (x_coord <= -0.6 || x_coord >= 0.6){
        alpha_o = 0.0;
    }
    else{
        alpha_o = 0.7;
    }
  }
  else{
    alpha_o = 0.0;
  }

  //initial shear modulus (which take initial damage into account)
  ADReal initial_shear_modulus = _shear_modulus_o + _xi_0 * alpha_o * _gamma_damaged_r;

  //compute initial strain based on initial stress
  /// lambda (first lame const)
  _lambda[_qp] = _lambda_o;
  /// mu (shear modulus)
  _shear_modulus[_qp] = initial_shear_modulus;
  /// gamma_damaged (damage modulus)
  _gamma_damaged[_qp] = alpha_o * _gamma_damaged_r;

  //allpha, B
  _alpha_damagedvar[_qp] = alpha_o;
  _B[_qp] = 0.0;

  //Convert (lambda_o,shear_modulus_o) to (youngs_modulus_o,poisson_ratio_o)
  ADReal youngs_modulus_o = initial_shear_modulus * ( 3 * _lambda_o + 2 * initial_shear_modulus ) / ( _lambda_o + initial_shear_modulus );
  ADReal poisson_ratio_o = _lambda_o / ( 2 * ( _lambda_o + initial_shear_modulus ));

  //Get stress components
  ADRankTwoTensor stress_initial = _static_initial_stress_tensor[_qp];

  //Get stress components
  ADReal sts11_init = stress_initial(0,0);
  ADReal sts12_init = stress_initial(0,1);
  ADReal sts22_init = stress_initial(1,1);
  //ADReal sts33_init = stress_initial(2,2);
  ADReal sts33_init = poisson_ratio_o * ( sts11_init + sts22_init );

  //Compute strain components using Hooke's Law
  ADReal eps11_init = 1.0 / youngs_modulus_o * ( sts11_init - poisson_ratio_o * ( sts22_init + sts33_init ) );
  ADReal eps22_init = 1.0 / youngs_modulus_o * ( sts22_init - poisson_ratio_o * ( sts11_init + sts33_init ) ); 
  ADReal eps12_init = 1.0 / youngs_modulus_o * ( ( 1 + poisson_ratio_o ) * sts12_init                       );
  ADReal eps13_init = 0.0;
  ADReal eps23_init = 0.0;
  ADReal eps33_init = 0.0;

  //Compute xi, I1, I2
  ADReal I1_init = eps11_init + eps22_init + eps33_init;
  ADReal I2_init = eps11_init * eps11_init + eps22_init * eps22_init + eps33_init * eps33_init + 2 * eps12_init * eps12_init + 2 * eps13_init * eps13_init + 2 * eps23_init * eps23_init;
  ADReal xi_init = I1_init / std::sqrt( I2_init );

  //Compute eps
  //eps_p
  _eps_p[_qp](0,0) = 0.0; _eps_p[_qp](0,1) = 0.0; _eps_p[_qp](0,2) = 0.0;
  _eps_p[_qp](1,0) = 0.0; _eps_p[_qp](1,1) = 0.0; _eps_p[_qp](1,2) = 0.0;
  _eps_p[_qp](2,0) = 0.0; _eps_p[_qp](2,1) = 0.0; _eps_p[_qp](2,2) = 0.0;
  //eps_e
  _eps_e[_qp](0,0) = eps11_init; _eps_e[_qp](0,1) = eps12_init; 
  _eps_e[_qp](1,0) = eps12_init; _eps_e[_qp](1,1) = eps22_init;
  _eps_e[_qp](2,2) = 0.0;
  //eps_total
  _eps_total[_qp](0,0) = eps11_init; _eps_total[_qp](0,1) = eps12_init;
  _eps_total[_qp](1,0) = eps12_init; _eps_total[_qp](1,1) = eps22_init;
  _eps_total[_qp](2,2) = 0.0;
  //sts_total
  _sts_total[_qp] = stress_initial;

  //I1
  _I1[_qp] = I1_init;
  //I2
  _I2[_qp] = I2_init;
  //xi
  _xi[_qp] = xi_init; 
}

//template class ADComputeDamageBreakageStressTempl<RankTwoTensor, RankFourTensor, ADReal>;
//template class ADComputeDamageBreakageStressTempl<SymmetricRankTwoTensor, SymmetricRankFourTensor>;