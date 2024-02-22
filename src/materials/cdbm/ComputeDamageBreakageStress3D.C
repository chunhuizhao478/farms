//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ComputeDamageBreakageStress3D.h"
#include "NestedSolve.h"
#include "FEProblem.h"

registerMooseObject("farmsApp", ComputeDamageBreakageStress3D);

InputParameters
ComputeDamageBreakageStress3D::validParams()
{ 
  //Note: lambda_o, shear_modulus_o is defined in "ComputeGeneralDamageBreakageStressBase"
  //to initialize _lambda, _shear_modulus material properties
  InputParameters params = ComputeDamageBreakageStressBase3D::validParams();
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
  params.addRequiredParam<Real>(              "m1", "coefficient of power law indexes");
  params.addRequiredParam<Real>(              "m2", "coefficient of power law indexes");
  params.addRequiredParam<Real>(               "D", "coefficient of alpha gradient");
  
  //variable parameters
  params.addRequiredCoupledVar("alpha_in", "damage variable computed from subApp");
  params.addRequiredCoupledVar(    "B_in", "breakage variable computed from subApp");
  params.addRequiredCoupledVar("alpha_grad_x", "damage variable gradient component in x computed from subApp");
  params.addRequiredCoupledVar("alpha_grad_y", "damage variable gradient component in y computed from subApp");
  params.addRequiredCoupledVar("alpha_grad_z", "damage variable gradient component in z computed from subApp");

  return params;
}

ComputeDamageBreakageStress3D::ComputeDamageBreakageStress3D(const InputParameters & parameters)
  : ComputeDamageBreakageStressBase3D(parameters),
    _static_initial_stress_tensor(getMaterialPropertyByName<RankTwoTensor>("static_initial_stress_tensor")),
    _static_initial_strain_tensor(getMaterialPropertyByName<RankTwoTensor>("static_initial_strain_tensor")),
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
    // _eqv_plastic_strain_old(getMaterialPropertyOldByName<Real>("eqv_plastic_strain")),
    _alpha_in(coupledValue("alpha_in")),
    _B_in(coupledValue("B_in")),
    _alpha_grad_x(coupledValue("alpha_grad_x")),
    _alpha_grad_y(coupledValue("alpha_grad_y")),
    _alpha_grad_z(coupledValue("alpha_grad_z")),
    _density_old(getMaterialPropertyOldByName<Real>("density")),
    _D(getParam<Real>("D"))
{
}

void
ComputeDamageBreakageStress3D::initialSetup()
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
ComputeDamageBreakageStress3D::computeQpStress()
{ 

  if (_t == 0.0)
  {
    setupInitial(); 
    //initialize lambda, shear_modulus, gamma_damaged
    //alpha_damagedvar, B, eps_p, eps_e, eps_total, sts_total, I1, I2, xi,
    //shear_wave_speed, pressure_wave_speed
  }
  else{

    //execute before system solve
    //if (_fe_problem.getCurrentExecuteOnFlag()=="LINEAR"){
      
      /* 
      compute alpha and B parameters
      */

      //alpha, B are updated
      Real alpha_out = _alpha_in[_qp];
      Real B_out = _B_in[_qp];

      //grad_alpha
      Real alpha_grad_x = _alpha_grad_x[_qp];
      Real alpha_grad_y = _alpha_grad_y[_qp];
      Real alpha_grad_z = _alpha_grad_z[_qp];
      Real D = _D;

      //save alpha and B
      _alpha_damagedvar[_qp] = alpha_out;
      _B[_qp] = B_out;

      /*
        update modulus
      */

      //lambda, shear_modulus, gamma_damaged are updated
      Real lambda_out = _lambda_o;
      Real shear_modulus_out = _shear_modulus_o + alpha_out * _xi_0 * _gamma_damaged_r;
      Real gamma_damaged_out = alpha_out * _gamma_damaged_r;

      //save
      _lambda[_qp] = lambda_out;
      _shear_modulus[_qp] = shear_modulus_out;
      _gamma_damaged[_qp] = gamma_damaged_out;

      /*
        update shear/pressure wave speed
      */

      // Real density = _density_old[_qp];

      // _shear_wave_speed[_qp] = std::sqrt( ( shear_modulus_out ) / ( density ) );
      // _pressure_wave_speed[_qp] = std::sqrt( ( lambda_out + 2 * shear_modulus_out ) / ( density ) );

      /*
        compute strain
      */

      //eps_p, eps_e, eps_total are updated
      //Retrieve parameter values
      //strain: total, viscoelastic
      RankTwoTensor eps_total_change = _mechanical_strain[_qp];
      RankTwoTensor eps_total_change_old = _mechanical_strain_old[_qp];
      RankTwoTensor eps_total_change_inc = eps_total_change - eps_total_change_old;

      //get total strain (add increment)
      RankTwoTensor eps_total = _eps_total_old[_qp] + eps_total_change_inc;

      //take the viscoelastic strain from the previous time step
      RankTwoTensor eps_viscoelastic = _eps_p_old[_qp];
        
      //take updated parameters
      Real lambda = lambda_out;
      Real shear_modulus = shear_modulus_out;
      Real gamma_damaged = gamma_damaged_out;
      Real B = B_out;

      //Define components
      Real eps11p_pre = eps_viscoelastic(0,0);
      Real eps22p_pre = eps_viscoelastic(1,1);
      Real eps33p_pre = eps_viscoelastic(2,2);
      Real eps12p_pre = eps_viscoelastic(0,1);
      Real eps13p_pre = eps_viscoelastic(0,2);
      Real eps23p_pre = eps_viscoelastic(1,2);
      Real eps11t = eps_total(0,0);
      Real eps22t = eps_total(1,1);
      Real eps33t = eps_total(2,2);
      Real eps12t = eps_total(0,1);
      Real eps13t = eps_total(0,2);
      Real eps23t = eps_total(1,2);

      auto compute = [&](const NestedSolve::Value<> & guess,
                          NestedSolve::Value<> & residual,
                          NestedSolve::Jacobian<> & jacobian)
      {
          //current viscoelastic strain (with increment unknown)
          //Note: 
          //guess(0) = eps11p_inc
          //guess(1) = eps22p_inc
          //guess(2) = eps33p_inc
          //guess(3) = eps12p_inc
          //guess(4) = eps13p_inc
          //guess(5) = eps23p_inc

          Real eps11p = eps11p_pre + guess(0);
          Real eps22p = eps22p_pre + guess(1);
          Real eps33p = eps33p_pre + guess(2);
          Real eps12p = eps12p_pre + guess(3);
          Real eps13p = eps13p_pre + guess(4);
          Real eps23p = eps23p_pre + guess(5);

          //compute elastic strain
          Real eps11e = eps11t - eps11p;
          Real eps22e = eps22t - eps22p;
          Real eps33e = eps33t - eps33p;
          Real eps12e = eps12t - eps12p;
          Real eps13e = eps13t - eps13p;
          Real eps23e = eps23t - eps23p;

          //represent I1 I2 xi using elastic strain
          Real I1 = eps11e + eps22e + eps33e;
          Real I2 = eps11e * eps11e + eps22e * eps22e + eps33e * eps33e + 2 * eps12e * eps12e + 2 * eps13e * eps13e + 2 * eps23e * eps23e;
          Real xi = I1 / sqrt(I2);

          //Represent sigma (solid(s) + granular(b))
          Real sigma11_s = ( lambda - gamma_damaged / xi ) * I1 + ( 2 * shear_modulus - gamma_damaged * xi ) * eps11e;
          Real sigma11_b = ( 2 * _a2 + _a1 / xi + 3 * _a3 * xi ) * I1 + ( 2 * _a0 + _a1 * xi - _a3 * pow(xi,3) ) * eps11e;
          Real sigma22_s = ( lambda - gamma_damaged / xi ) * I1 + ( 2 * shear_modulus - gamma_damaged * xi ) * eps22e;
          Real sigma22_b = ( 2 * _a2 + _a1 / xi + 3 * _a3 * xi ) * I1 + ( 2 * _a0 + _a1 * xi - _a3 * pow(xi,3) ) * eps22e;
          Real sigma33_s = ( lambda - gamma_damaged / xi ) * I1 + ( 2 * shear_modulus - gamma_damaged * xi ) * eps33e;
          Real sigma33_b = ( 2 * _a2 + _a1 / xi + 3 * _a3 * xi ) * I1 + ( 2 * _a0 + _a1 * xi - _a3 * pow(xi,3) ) * eps33e;
          Real sigma12_s = ( 2 * shear_modulus - gamma_damaged * xi ) * eps12e;
          Real sigma12_b = ( 2 * _a0 + _a1 * xi - _a3 * pow(xi,3) ) * eps12e;
          Real sigma13_s = ( 2 * shear_modulus - gamma_damaged * xi ) * eps13e;
          Real sigma13_b = ( 2 * _a0 + _a1 * xi - _a3 * pow(xi,3) ) * eps13e;
          Real sigma23_s = ( 2 * shear_modulus - gamma_damaged * xi ) * eps23e;
          Real sigma23_b = ( 2 * _a0 + _a1 * xi - _a3 * pow(xi,3) ) * eps23e;

          //Represent total stress
          Real sigma11_t = (1 - B) * sigma11_s + B * sigma11_b;
          Real sigma22_t = (1 - B) * sigma22_s + B * sigma22_b;
          Real sigma33_t = (1 - B) * sigma33_s + B * sigma33_b;
          Real sigma12_t = (1 - B) * sigma12_s + B * sigma12_b;
          Real sigma13_t = (1 - B) * sigma13_s + B * sigma13_b;
          Real sigma23_t = (1 - B) * sigma23_s + B * sigma23_b;

          //Represent deviatoric stress
          Real sigma_d11 = sigma11_t - 1/3 * (sigma11_t + sigma22_t + sigma33_t);
          Real sigma_d22 = sigma22_t - 1/3 * (sigma11_t + sigma22_t + sigma33_t);
          Real sigma_d33 = sigma33_t - 1/3 * (sigma11_t + sigma22_t + sigma33_t);
          Real sigma_d12 = sigma12_t;
          Real sigma_d13 = sigma13_t;
          Real sigma_d23 = sigma23_t;

          //Setup equations
          ///residual
          residual(0) = guess(0) - _dt * _C_g * pow(B,_m1) * pow(sigma_d11,_m2);
          residual(1) = guess(1) - _dt * _C_g * pow(B,_m1) * pow(sigma_d22,_m2);
          residual(2) = guess(2) - _dt * _C_g * pow(B,_m1) * pow(sigma_d33,_m2);
          residual(3) = guess(3) - _dt * _C_g * pow(B,_m1) * pow(sigma_d12,_m2);
          residual(4) = guess(4) - _dt * _C_g * pow(B,_m1) * pow(sigma_d13,_m2);
          residual(5) = guess(5) - _dt * _C_g * pow(B,_m1) * pow(sigma_d23,_m2);

          ///jacobian //too long ...
          ///Define parameters
          Real minus_eps11e = guess(0) - eps11t + eps11p_pre;
          Real minus_eps22e = guess(1) - eps22t + eps22p_pre;
          Real minus_eps33e = guess(2) - eps33t + eps33p_pre;
          Real minus_eps12e = guess(3) - eps12t + eps12p_pre;
          Real minus_eps13e = guess(4) - eps13t + eps13p_pre;
          Real minus_eps23e = guess(5) - eps23t + eps23p_pre;
          Real minus_eps_tr = minus_eps11e + minus_eps22e + minus_eps33e;
          Real minus_epse_I2 = pow(minus_eps11e,2) + 2*pow(minus_eps12e,2) + pow(minus_eps22e,2) + pow(minus_eps33e,2) + 2*pow(minus_eps13e,2) + 2*pow(minus_eps23e,2);
          jacobian(0,0) = pow(B,_m1)*_C_g*_dt*_m2*pow((2*((2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5))*(minus_eps11e) + (lambda + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr))*(B - 1))/3 - (((2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5))*(minus_eps22e) + (lambda + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr))*(B - 1))/3 - (((2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5))*(minus_eps33e) + (lambda + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr))*(B - 1))/3 - (2*B*((2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5))*(minus_eps11e) - ((3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) - 2*_a2 + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr)))/3 + (B*((2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5))*(minus_eps22e) - ((3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) - 2*_a2 + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr)))/3 + (B*((2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5))*(minus_eps33e) - ((3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) - 2*_a2 + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr)))/3,_m2 - 1)*(((B - 1)*(lambda - ((gamma_damaged*pow(minus_epse_I2,0.5))/pow(minus_eps_tr,2.0) - (gamma_damaged*(2*minus_eps11e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr) + (gamma_damaged/pow(minus_epse_I2,0.5) - (gamma_damaged*(2*minus_eps11e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)))*(minus_eps22e) + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr)))/3 - (2*B*((minus_eps11e)*(_a1/pow(minus_epse_I2,0.5) - (3*_a3*pow(minus_eps_tr,2.0))/pow(minus_epse_I2,1.5) - (_a1*(2*minus_eps11e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) + (3*_a3*(2*minus_eps11e)*pow(minus_eps_tr,3.0))/(2*pow(minus_epse_I2,2.5))) - 2*_a2 - 2*_a0 + ((3*_a3)/pow(minus_epse_I2,0.5) - (_a1*pow(minus_epse_I2,0.5))/pow(minus_eps_tr,2.0) - (3*_a3*(2*minus_eps11e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) + (_a1*(2*minus_eps11e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr) + (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr) - (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5)))/3 + ((B - 1)*(lambda - ((gamma_damaged*pow(minus_epse_I2,0.5))/pow(minus_eps_tr,2.0) - (gamma_damaged*(2*minus_eps11e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr) + (gamma_damaged/pow(minus_epse_I2,0.5) - (gamma_damaged*(2*minus_eps11e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)))*(minus_eps33e) + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr)))/3 - (2*(B - 1)*(lambda + 2*shear_modulus - ((gamma_damaged*pow(minus_epse_I2,0.5))/pow(minus_eps_tr,2.0) - (gamma_damaged*(2*minus_eps11e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr) + (gamma_damaged/pow(minus_epse_I2,0.5) - (gamma_damaged*(2*minus_eps11e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)))*(minus_eps11e) + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr)))/3 + (B*((minus_eps22e)*(_a1/pow(minus_epse_I2,0.5) - (3*_a3*pow(minus_eps_tr,2.0))/pow(minus_epse_I2,1.5) - (_a1*(2*minus_eps11e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) + (3*_a3*(2*minus_eps11e)*pow(minus_eps_tr,3.0))/(2*pow(minus_epse_I2,2.5))) - 2*_a2 + ((3*_a3)/pow(minus_epse_I2,0.5) - (_a1*pow(minus_epse_I2,0.5))/pow(minus_eps_tr,2.0) - (3*_a3*(2*minus_eps11e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) + (_a1*(2*minus_eps11e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr) + (3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr)))/3 + (B*((minus_eps33e)*(_a1/pow(minus_epse_I2,0.5) - (3*_a3*pow(minus_eps_tr,2.0))/pow(minus_epse_I2,1.5) - (_a1*(2*minus_eps11e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) + (3*_a3*(2*minus_eps11e)*pow(minus_eps_tr,3.0))/(2*pow(minus_epse_I2,2.5))) - 2*_a2 + ((3*_a3)/pow(minus_epse_I2,0.5) - (_a1*pow(minus_epse_I2,0.5))/pow(minus_eps_tr,2.0) - (3*_a3*(2*minus_eps11e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) + (_a1*(2*minus_eps11e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr) + (3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr)))/3) + 1;
          jacobian(0,1) = pow(B,_m1)*_C_g*_dt*_m2*pow((2*((2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5))*(minus_eps11e) + (lambda + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr))*(B - 1))/3 - (((2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5))*(minus_eps22e) + (lambda + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr))*(B - 1))/3 - (((2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5))*(minus_eps33e) + (lambda + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr))*(B - 1))/3 - (2*B*((2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5))*(minus_eps11e) - ((3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) - 2*_a2 + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr)))/3 + (B*((2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5))*(minus_eps22e) - ((3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) - 2*_a2 + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr)))/3 + (B*((2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5))*(minus_eps33e) - ((3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) - 2*_a2 + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr)))/3,_m2 - 1)*((B*((minus_eps22e)*(_a1/pow(minus_epse_I2,0.5) - (3*_a3*pow(minus_eps_tr,2.0))/pow(minus_epse_I2,1.5) - (_a1*(2*minus_eps22e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) + (3*_a3*(2*minus_eps22e)*pow(minus_eps_tr,3.0))/(2*pow(minus_epse_I2,2.5))) - 2*_a2 - 2*_a0 + ((3*_a3)/pow(minus_epse_I2,0.5) - (_a1*pow(minus_epse_I2,0.5))/pow(minus_eps_tr,2.0) - (3*_a3*(2*minus_eps22e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) + (_a1*(2*minus_eps22e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr) + (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr) - (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5)))/3 - (2*(B - 1)*(lambda - ((gamma_damaged*pow(minus_epse_I2,0.5))/pow(minus_eps_tr,2.0) - (gamma_damaged*(2*minus_eps22e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr) + (gamma_damaged/pow(minus_epse_I2,0.5) - (gamma_damaged*(2*minus_eps22e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)))*(minus_eps11e) + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr)))/3 + ((B - 1)*(lambda - ((gamma_damaged*pow(minus_epse_I2,0.5))/pow(minus_eps_tr,2.0) - (gamma_damaged*(2*minus_eps22e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr) + (gamma_damaged/pow(minus_epse_I2,0.5) - (gamma_damaged*(2*minus_eps22e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)))*(minus_eps33e) + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr)))/3 + ((B - 1)*(lambda + 2*shear_modulus - ((gamma_damaged*pow(minus_epse_I2,0.5))/pow(minus_eps_tr,2.0) - (gamma_damaged*(2*minus_eps22e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr) + (gamma_damaged/pow(minus_epse_I2,0.5) - (gamma_damaged*(2*minus_eps22e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)))*(minus_eps22e) + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr)))/3 - (2*B*((minus_eps11e)*(_a1/pow(minus_epse_I2,0.5) - (3*_a3*pow(minus_eps_tr,2.0))/pow(minus_epse_I2,1.5) - (_a1*(2*minus_eps22e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) + (3*_a3*(2*minus_eps22e)*pow(minus_eps_tr,3.0))/(2*pow(minus_epse_I2,2.5))) - 2*_a2 + ((3*_a3)/pow(minus_epse_I2,0.5) - (_a1*pow(minus_epse_I2,0.5))/pow(minus_eps_tr,2.0) - (3*_a3*(2*minus_eps22e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) + (_a1*(2*minus_eps22e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr) + (3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr)))/3 + (B*((minus_eps33e)*(_a1/pow(minus_epse_I2,0.5) - (3*_a3*pow(minus_eps_tr,2.0))/pow(minus_epse_I2,1.5) - (_a1*(2*minus_eps22e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) + (3*_a3*(2*minus_eps22e)*pow(minus_eps_tr,3.0))/(2*pow(minus_epse_I2,2.5))) - 2*_a2 + ((3*_a3)/pow(minus_epse_I2,0.5) - (_a1*pow(minus_epse_I2,0.5))/pow(minus_eps_tr,2.0) - (3*_a3*(2*minus_eps22e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) + (_a1*(2*minus_eps22e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr) + (3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr)))/3);
          jacobian(0,2) = pow(B,_m1)*_C_g*_dt*_m2*pow((2*((2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5))*(minus_eps11e) + (lambda + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr))*(B - 1))/3 - (((2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5))*(minus_eps22e) + (lambda + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr))*(B - 1))/3 - (((2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5))*(minus_eps33e) + (lambda + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr))*(B - 1))/3 - (2*B*((2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5))*(minus_eps11e) - ((3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) - 2*_a2 + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr)))/3 + (B*((2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5))*(minus_eps22e) - ((3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) - 2*_a2 + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr)))/3 + (B*((2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5))*(minus_eps33e) - ((3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) - 2*_a2 + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr)))/3,_m2 - 1)*((B*((minus_eps33e)*(_a1/pow(minus_epse_I2,0.5) - (3*_a3*pow(minus_eps_tr,2.0))/pow(minus_epse_I2,1.5) - (_a1*(2*minus_eps33e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) + (3*_a3*(2*minus_eps33e)*pow(minus_eps_tr,3.0))/(2*pow(minus_epse_I2,2.5))) - 2*_a2 - 2*_a0 + ((3*_a3)/pow(minus_epse_I2,0.5) - (_a1*pow(minus_epse_I2,0.5))/pow(minus_eps_tr,2.0) - (3*_a3*(2*minus_eps33e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) + (_a1*(2*minus_eps33e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr) + (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr) - (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5)))/3 - (2*(B - 1)*(lambda - ((gamma_damaged*pow(minus_epse_I2,0.5))/pow(minus_eps_tr,2.0) - (gamma_damaged*(2*minus_eps33e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr) + (gamma_damaged/pow(minus_epse_I2,0.5) - (gamma_damaged*(2*minus_eps33e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)))*(minus_eps11e) + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr)))/3 + ((B - 1)*(lambda - ((gamma_damaged*pow(minus_epse_I2,0.5))/pow(minus_eps_tr,2.0) - (gamma_damaged*(2*minus_eps33e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr) + (gamma_damaged/pow(minus_epse_I2,0.5) - (gamma_damaged*(2*minus_eps33e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)))*(minus_eps22e) + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr)))/3 + ((B - 1)*(lambda + 2*shear_modulus - ((gamma_damaged*pow(minus_epse_I2,0.5))/pow(minus_eps_tr,2.0) - (gamma_damaged*(2*minus_eps33e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr) + (gamma_damaged/pow(minus_epse_I2,0.5) - (gamma_damaged*(2*minus_eps33e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)))*(minus_eps33e) + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr)))/3 - (2*B*((minus_eps11e)*(_a1/pow(minus_epse_I2,0.5) - (3*_a3*pow(minus_eps_tr,2.0))/pow(minus_epse_I2,1.5) - (_a1*(2*minus_eps33e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) + (3*_a3*(2*minus_eps33e)*pow(minus_eps_tr,3.0))/(2*pow(minus_epse_I2,2.5))) - 2*_a2 + ((3*_a3)/pow(minus_epse_I2,0.5) - (_a1*pow(minus_epse_I2,0.5))/pow(minus_eps_tr,2.0) - (3*_a3*(2*minus_eps33e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) + (_a1*(2*minus_eps33e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr) + (3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr)))/3 + (B*((minus_eps22e)*(_a1/pow(minus_epse_I2,0.5) - (3*_a3*pow(minus_eps_tr,2.0))/pow(minus_epse_I2,1.5) - (_a1*(2*minus_eps33e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) + (3*_a3*(2*minus_eps33e)*pow(minus_eps_tr,3.0))/(2*pow(minus_epse_I2,2.5))) - 2*_a2 + ((3*_a3)/pow(minus_epse_I2,0.5) - (_a1*pow(minus_epse_I2,0.5))/pow(minus_eps_tr,2.0) - (3*_a3*(2*minus_eps33e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) + (_a1*(2*minus_eps33e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr) + (3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr)))/3);
          jacobian(0,3) = pow(B,_m1)*_C_g*_dt*_m2*pow((2*((2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5))*(minus_eps11e) + (lambda + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr))*(B - 1))/3 - (((2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5))*(minus_eps22e) + (lambda + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr))*(B - 1))/3 - (((2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5))*(minus_eps33e) + (lambda + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr))*(B - 1))/3 - (2*B*((2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5))*(minus_eps11e) - ((3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) - 2*_a2 + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr)))/3 + (B*((2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5))*(minus_eps22e) - ((3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) - 2*_a2 + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr)))/3 + (B*((2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5))*(minus_eps33e) - ((3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) - 2*_a2 + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr)))/3,_m2 - 1)*((2*B*(((_a1*(4*minus_eps12e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) - (3*_a3*(4*minus_eps12e)*pow(minus_eps_tr,3.0))/(2*pow(minus_epse_I2,2.5)))*(minus_eps11e) + ((3*_a3*(4*minus_eps12e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) - (_a1*(4*minus_eps12e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr)))/3 - (B*(((_a1*(4*minus_eps12e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) - (3*_a3*(4*minus_eps12e)*pow(minus_eps_tr,3.0))/(2*pow(minus_epse_I2,2.5)))*(minus_eps22e) + ((3*_a3*(4*minus_eps12e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) - (_a1*(4*minus_eps12e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr)))/3 - (B*(((_a1*(4*minus_eps12e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) - (3*_a3*(4*minus_eps12e)*pow(minus_eps_tr,3.0))/(2*pow(minus_epse_I2,2.5)))*(minus_eps33e) + ((3*_a3*(4*minus_eps12e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) - (_a1*(4*minus_eps12e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr)))/3 - (2*(B - 1)*((gamma_damaged*(4*minus_eps12e))/(2*pow(minus_epse_I2,0.5)) - (gamma_damaged*(4*minus_eps12e)*(minus_eps11e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5))))/3 + ((B - 1)*((gamma_damaged*(4*minus_eps12e))/(2*pow(minus_epse_I2,0.5)) - (gamma_damaged*(4*minus_eps12e)*(minus_eps22e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5))))/3 + ((B - 1)*((gamma_damaged*(4*minus_eps12e))/(2*pow(minus_epse_I2,0.5)) - (gamma_damaged*(4*minus_eps12e)*(minus_eps33e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5))))/3);
          jacobian(0,4) = pow(B,_m1)*_C_g*_dt*_m2*pow((2*((2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5))*(minus_eps11e) + (lambda + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr))*(B - 1))/3 - (((2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5))*(minus_eps22e) + (lambda + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr))*(B - 1))/3 - (((2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5))*(minus_eps33e) + (lambda + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr))*(B - 1))/3 - (2*B*((2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5))*(minus_eps11e) - ((3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) - 2*_a2 + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr)))/3 + (B*((2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5))*(minus_eps22e) - ((3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) - 2*_a2 + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr)))/3 + (B*((2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5))*(minus_eps33e) - ((3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) - 2*_a2 + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr)))/3,_m2 - 1)*((2*B*(((_a1*(4*minus_eps13e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) - (3*_a3*(4*minus_eps13e)*pow(minus_eps_tr,3.0))/(2*pow(minus_epse_I2,2.5)))*(minus_eps11e) + ((3*_a3*(4*minus_eps13e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) - (_a1*(4*minus_eps13e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr)))/3 - (B*(((_a1*(4*minus_eps13e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) - (3*_a3*(4*minus_eps13e)*pow(minus_eps_tr,3.0))/(2*pow(minus_epse_I2,2.5)))*(minus_eps22e) + ((3*_a3*(4*minus_eps13e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) - (_a1*(4*minus_eps13e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr)))/3 - (B*(((_a1*(4*minus_eps13e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) - (3*_a3*(4*minus_eps13e)*pow(minus_eps_tr,3.0))/(2*pow(minus_epse_I2,2.5)))*(minus_eps33e) + ((3*_a3*(4*minus_eps13e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) - (_a1*(4*minus_eps13e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr)))/3 - (2*(B - 1)*((gamma_damaged*(4*minus_eps13e))/(2*pow(minus_epse_I2,0.5)) - (gamma_damaged*(4*minus_eps13e)*(minus_eps11e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5))))/3 + ((B - 1)*((gamma_damaged*(4*minus_eps13e))/(2*pow(minus_epse_I2,0.5)) - (gamma_damaged*(4*minus_eps13e)*(minus_eps22e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5))))/3 + ((B - 1)*((gamma_damaged*(4*minus_eps13e))/(2*pow(minus_epse_I2,0.5)) - (gamma_damaged*(4*minus_eps13e)*(minus_eps33e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5))))/3);
          jacobian(0,5) = pow(B,_m1)*_C_g*_dt*_m2*pow((2*((2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5))*(minus_eps11e) + (lambda + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr))*(B - 1))/3 - (((2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5))*(minus_eps22e) + (lambda + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr))*(B - 1))/3 - (((2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5))*(minus_eps33e) + (lambda + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr))*(B - 1))/3 - (2*B*((2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5))*(minus_eps11e) - ((3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) - 2*_a2 + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr)))/3 + (B*((2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5))*(minus_eps22e) - ((3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) - 2*_a2 + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr)))/3 + (B*((2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5))*(minus_eps33e) - ((3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) - 2*_a2 + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr)))/3,_m2 - 1)*((2*B*(((_a1*(4*minus_eps23e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) - (3*_a3*(4*minus_eps23e)*pow(minus_eps_tr,3.0))/(2*pow(minus_epse_I2,2.5)))*(minus_eps11e) + ((3*_a3*(4*minus_eps23e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) - (_a1*(4*minus_eps23e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr)))/3 - (B*(((_a1*(4*minus_eps23e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) - (3*_a3*(4*minus_eps23e)*pow(minus_eps_tr,3.0))/(2*pow(minus_epse_I2,2.5)))*(minus_eps22e) + ((3*_a3*(4*minus_eps23e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) - (_a1*(4*minus_eps23e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr)))/3 - (B*(((_a1*(4*minus_eps23e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) - (3*_a3*(4*minus_eps23e)*pow(minus_eps_tr,3.0))/(2*pow(minus_epse_I2,2.5)))*(minus_eps33e) + ((3*_a3*(4*minus_eps23e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) - (_a1*(4*minus_eps23e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr)))/3 - (2*(B - 1)*((gamma_damaged*(4*minus_eps23e))/(2*pow(minus_epse_I2,0.5)) - (gamma_damaged*(4*minus_eps23e)*(minus_eps11e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5))))/3 + ((B - 1)*((gamma_damaged*(4*minus_eps23e))/(2*pow(minus_epse_I2,0.5)) - (gamma_damaged*(4*minus_eps23e)*(minus_eps22e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5))))/3 + ((B - 1)*((gamma_damaged*(4*minus_eps23e))/(2*pow(minus_epse_I2,0.5)) - (gamma_damaged*(4*minus_eps23e)*(minus_eps33e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5))))/3);
          jacobian(1,0) = pow(B,_m1)*_C_g*_dt*_m2*pow((2*((2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5))*(minus_eps22e) + (lambda + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr))*(B - 1))/3 - (((2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5))*(minus_eps11e) + (lambda + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr))*(B - 1))/3 - (((2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5))*(minus_eps33e) + (lambda + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr))*(B - 1))/3 + (B*((2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5))*(minus_eps11e) - ((3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) - 2*_a2 + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr)))/3 - (2*B*((2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5))*(minus_eps22e) - ((3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) - 2*_a2 + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr)))/3 + (B*((2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5))*(minus_eps33e) - ((3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) - 2*_a2 + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr)))/3,_m2 - 1)*((B*((minus_eps11e)*(_a1/pow(minus_epse_I2,0.5) - (3*_a3*pow(minus_eps_tr,2.0))/pow(minus_epse_I2,1.5) - (_a1*(2*minus_eps11e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) + (3*_a3*(2*minus_eps11e)*pow(minus_eps_tr,3.0))/(2*pow(minus_epse_I2,2.5))) - 2*_a2 - 2*_a0 + ((3*_a3)/pow(minus_epse_I2,0.5) - (_a1*pow(minus_epse_I2,0.5))/pow(minus_eps_tr,2.0) - (3*_a3*(2*minus_eps11e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) + (_a1*(2*minus_eps11e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr) + (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr) - (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5)))/3 - (2*(B - 1)*(lambda - ((gamma_damaged*pow(minus_epse_I2,0.5))/pow(minus_eps_tr,2.0) - (gamma_damaged*(2*minus_eps11e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr) + (gamma_damaged/pow(minus_epse_I2,0.5) - (gamma_damaged*(2*minus_eps11e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)))*(minus_eps22e) + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr)))/3 + ((B - 1)*(lambda - ((gamma_damaged*pow(minus_epse_I2,0.5))/pow(minus_eps_tr,2.0) - (gamma_damaged*(2*minus_eps11e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr) + (gamma_damaged/pow(minus_epse_I2,0.5) - (gamma_damaged*(2*minus_eps11e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)))*(minus_eps33e) + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr)))/3 + ((B - 1)*(lambda + 2*shear_modulus - ((gamma_damaged*pow(minus_epse_I2,0.5))/pow(minus_eps_tr,2.0) - (gamma_damaged*(2*minus_eps11e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr) + (gamma_damaged/pow(minus_epse_I2,0.5) - (gamma_damaged*(2*minus_eps11e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)))*(minus_eps11e) + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr)))/3 - (2*B*((minus_eps22e)*(_a1/pow(minus_epse_I2,0.5) - (3*_a3*pow(minus_eps_tr,2.0))/pow(minus_epse_I2,1.5) - (_a1*(2*minus_eps11e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) + (3*_a3*(2*minus_eps11e)*pow(minus_eps_tr,3.0))/(2*pow(minus_epse_I2,2.5))) - 2*_a2 + ((3*_a3)/pow(minus_epse_I2,0.5) - (_a1*pow(minus_epse_I2,0.5))/pow(minus_eps_tr,2.0) - (3*_a3*(2*minus_eps11e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) + (_a1*(2*minus_eps11e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr) + (3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr)))/3 + (B*((minus_eps33e)*(_a1/pow(minus_epse_I2,0.5) - (3*_a3*pow(minus_eps_tr,2.0))/pow(minus_epse_I2,1.5) - (_a1*(2*minus_eps11e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) + (3*_a3*(2*minus_eps11e)*pow(minus_eps_tr,3.0))/(2*pow(minus_epse_I2,2.5))) - 2*_a2 + ((3*_a3)/pow(minus_epse_I2,0.5) - (_a1*pow(minus_epse_I2,0.5))/pow(minus_eps_tr,2.0) - (3*_a3*(2*minus_eps11e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) + (_a1*(2*minus_eps11e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr) + (3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr)))/3);
          jacobian(1,1) = pow(B,_m1)*_C_g*_dt*_m2*pow((2*((2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5))*(minus_eps22e) + (lambda + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr))*(B - 1))/3 - (((2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5))*(minus_eps11e) + (lambda + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr))*(B - 1))/3 - (((2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5))*(minus_eps33e) + (lambda + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr))*(B - 1))/3 + (B*((2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5))*(minus_eps11e) - ((3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) - 2*_a2 + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr)))/3 - (2*B*((2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5))*(minus_eps22e) - ((3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) - 2*_a2 + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr)))/3 + (B*((2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5))*(minus_eps33e) - ((3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) - 2*_a2 + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr)))/3,_m2 - 1)*(((B - 1)*(lambda - ((gamma_damaged*pow(minus_epse_I2,0.5))/pow(minus_eps_tr,2.0) - (gamma_damaged*(2*minus_eps22e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr) + (gamma_damaged/pow(minus_epse_I2,0.5) - (gamma_damaged*(2*minus_eps22e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)))*(minus_eps11e) + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr)))/3 - (2*B*((minus_eps22e)*(_a1/pow(minus_epse_I2,0.5) - (3*_a3*pow(minus_eps_tr,2.0))/pow(minus_epse_I2,1.5) - (_a1*(2*minus_eps22e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) + (3*_a3*(2*minus_eps22e)*pow(minus_eps_tr,3.0))/(2*pow(minus_epse_I2,2.5))) - 2*_a2 - 2*_a0 + ((3*_a3)/pow(minus_epse_I2,0.5) - (_a1*pow(minus_epse_I2,0.5))/pow(minus_eps_tr,2.0) - (3*_a3*(2*minus_eps22e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) + (_a1*(2*minus_eps22e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr) + (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr) - (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5)))/3 + ((B - 1)*(lambda - ((gamma_damaged*pow(minus_epse_I2,0.5))/pow(minus_eps_tr,2.0) - (gamma_damaged*(2*minus_eps22e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr) + (gamma_damaged/pow(minus_epse_I2,0.5) - (gamma_damaged*(2*minus_eps22e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)))*(minus_eps33e) + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr)))/3 - (2*(B - 1)*(lambda + 2*shear_modulus - ((gamma_damaged*pow(minus_epse_I2,0.5))/pow(minus_eps_tr,2.0) - (gamma_damaged*(2*minus_eps22e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr) + (gamma_damaged/pow(minus_epse_I2,0.5) - (gamma_damaged*(2*minus_eps22e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)))*(minus_eps22e) + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr)))/3 + (B*((minus_eps11e)*(_a1/pow(minus_epse_I2,0.5) - (3*_a3*pow(minus_eps_tr,2.0))/pow(minus_epse_I2,1.5) - (_a1*(2*minus_eps22e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) + (3*_a3*(2*minus_eps22e)*pow(minus_eps_tr,3.0))/(2*pow(minus_epse_I2,2.5))) - 2*_a2 + ((3*_a3)/pow(minus_epse_I2,0.5) - (_a1*pow(minus_epse_I2,0.5))/pow(minus_eps_tr,2.0) - (3*_a3*(2*minus_eps22e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) + (_a1*(2*minus_eps22e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr) + (3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr)))/3 + (B*((minus_eps33e)*(_a1/pow(minus_epse_I2,0.5) - (3*_a3*pow(minus_eps_tr,2.0))/pow(minus_epse_I2,1.5) - (_a1*(2*minus_eps22e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) + (3*_a3*(2*minus_eps22e)*pow(minus_eps_tr,3.0))/(2*pow(minus_epse_I2,2.5))) - 2*_a2 + ((3*_a3)/pow(minus_epse_I2,0.5) - (_a1*pow(minus_epse_I2,0.5))/pow(minus_eps_tr,2.0) - (3*_a3*(2*minus_eps22e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) + (_a1*(2*minus_eps22e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr) + (3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr)))/3) + 1;
          jacobian(1,2) = pow(B,_m1)*_C_g*_dt*_m2*pow((2*((2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5))*(minus_eps22e) + (lambda + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr))*(B - 1))/3 - (((2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5))*(minus_eps11e) + (lambda + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr))*(B - 1))/3 - (((2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5))*(minus_eps33e) + (lambda + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr))*(B - 1))/3 + (B*((2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5))*(minus_eps11e) - ((3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) - 2*_a2 + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr)))/3 - (2*B*((2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5))*(minus_eps22e) - ((3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) - 2*_a2 + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr)))/3 + (B*((2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5))*(minus_eps33e) - ((3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) - 2*_a2 + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr)))/3,_m2 - 1)*((B*((minus_eps33e)*(_a1/pow(minus_epse_I2,0.5) - (3*_a3*pow(minus_eps_tr,2.0))/pow(minus_epse_I2,1.5) - (_a1*(2*minus_eps33e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) + (3*_a3*(2*minus_eps33e)*pow(minus_eps_tr,3.0))/(2*pow(minus_epse_I2,2.5))) - 2*_a2 - 2*_a0 + ((3*_a3)/pow(minus_epse_I2,0.5) - (_a1*pow(minus_epse_I2,0.5))/pow(minus_eps_tr,2.0) - (3*_a3*(2*minus_eps33e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) + (_a1*(2*minus_eps33e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr) + (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr) - (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5)))/3 + ((B - 1)*(lambda - ((gamma_damaged*pow(minus_epse_I2,0.5))/pow(minus_eps_tr,2.0) - (gamma_damaged*(2*minus_eps33e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr) + (gamma_damaged/pow(minus_epse_I2,0.5) - (gamma_damaged*(2*minus_eps33e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)))*(minus_eps11e) + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr)))/3 - (2*(B - 1)*(lambda - ((gamma_damaged*pow(minus_epse_I2,0.5))/pow(minus_eps_tr,2.0) - (gamma_damaged*(2*minus_eps33e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr) + (gamma_damaged/pow(minus_epse_I2,0.5) - (gamma_damaged*(2*minus_eps33e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)))*(minus_eps22e) + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr)))/3 + ((B - 1)*(lambda + 2*shear_modulus - ((gamma_damaged*pow(minus_epse_I2,0.5))/pow(minus_eps_tr,2.0) - (gamma_damaged*(2*minus_eps33e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr) + (gamma_damaged/pow(minus_epse_I2,0.5) - (gamma_damaged*(2*minus_eps33e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)))*(minus_eps33e) + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr)))/3 + (B*((minus_eps11e)*(_a1/pow(minus_epse_I2,0.5) - (3*_a3*pow(minus_eps_tr,2.0))/pow(minus_epse_I2,1.5) - (_a1*(2*minus_eps33e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) + (3*_a3*(2*minus_eps33e)*pow(minus_eps_tr,3.0))/(2*pow(minus_epse_I2,2.5))) - 2*_a2 + ((3*_a3)/pow(minus_epse_I2,0.5) - (_a1*pow(minus_epse_I2,0.5))/pow(minus_eps_tr,2.0) - (3*_a3*(2*minus_eps33e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) + (_a1*(2*minus_eps33e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr) + (3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr)))/3 - (2*B*((minus_eps22e)*(_a1/pow(minus_epse_I2,0.5) - (3*_a3*pow(minus_eps_tr,2.0))/pow(minus_epse_I2,1.5) - (_a1*(2*minus_eps33e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) + (3*_a3*(2*minus_eps33e)*pow(minus_eps_tr,3.0))/(2*pow(minus_epse_I2,2.5))) - 2*_a2 + ((3*_a3)/pow(minus_epse_I2,0.5) - (_a1*pow(minus_epse_I2,0.5))/pow(minus_eps_tr,2.0) - (3*_a3*(2*minus_eps33e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) + (_a1*(2*minus_eps33e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr) + (3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr)))/3);
          jacobian(1,3) = -pow(B,_m1)*_C_g*_dt*_m2*pow((2*((2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5))*(minus_eps22e) + (lambda + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr))*(B - 1))/3 - (((2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5))*(minus_eps11e) + (lambda + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr))*(B - 1))/3 - (((2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5))*(minus_eps33e) + (lambda + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr))*(B - 1))/3 + (B*((2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5))*(minus_eps11e) - ((3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) - 2*_a2 + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr)))/3 - (2*B*((2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5))*(minus_eps22e) - ((3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) - 2*_a2 + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr)))/3 + (B*((2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5))*(minus_eps33e) - ((3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) - 2*_a2 + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr)))/3,_m2 - 1)*((B*(((_a1*(4*minus_eps12e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) - (3*_a3*(4*minus_eps12e)*pow(minus_eps_tr,3.0))/(2*pow(minus_epse_I2,2.5)))*(minus_eps11e) + ((3*_a3*(4*minus_eps12e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) - (_a1*(4*minus_eps12e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr)))/3 - (2*B*(((_a1*(4*minus_eps12e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) - (3*_a3*(4*minus_eps12e)*pow(minus_eps_tr,3.0))/(2*pow(minus_epse_I2,2.5)))*(minus_eps22e) + ((3*_a3*(4*minus_eps12e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) - (_a1*(4*minus_eps12e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr)))/3 + (B*(((_a1*(4*minus_eps12e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) - (3*_a3*(4*minus_eps12e)*pow(minus_eps_tr,3.0))/(2*pow(minus_epse_I2,2.5)))*(minus_eps33e) + ((3*_a3*(4*minus_eps12e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) - (_a1*(4*minus_eps12e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr)))/3 - ((B - 1)*((gamma_damaged*(4*minus_eps12e))/(2*pow(minus_epse_I2,0.5)) - (gamma_damaged*(4*minus_eps12e)*(minus_eps11e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5))))/3 + (2*(B - 1)*((gamma_damaged*(4*minus_eps12e))/(2*pow(minus_epse_I2,0.5)) - (gamma_damaged*(4*minus_eps12e)*(minus_eps22e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5))))/3 - ((B - 1)*((gamma_damaged*(4*minus_eps12e))/(2*pow(minus_epse_I2,0.5)) - (gamma_damaged*(4*minus_eps12e)*(minus_eps33e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5))))/3);
          jacobian(1,4) = -pow(B,_m1)*_C_g*_dt*_m2*pow((2*((2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5))*(minus_eps22e) + (lambda + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr))*(B - 1))/3 - (((2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5))*(minus_eps11e) + (lambda + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr))*(B - 1))/3 - (((2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5))*(minus_eps33e) + (lambda + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr))*(B - 1))/3 + (B*((2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5))*(minus_eps11e) - ((3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) - 2*_a2 + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr)))/3 - (2*B*((2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5))*(minus_eps22e) - ((3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) - 2*_a2 + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr)))/3 + (B*((2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5))*(minus_eps33e) - ((3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) - 2*_a2 + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr)))/3,_m2 - 1)*((B*(((_a1*(4*minus_eps13e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) - (3*_a3*(4*minus_eps13e)*pow(minus_eps_tr,3.0))/(2*pow(minus_epse_I2,2.5)))*(minus_eps11e) + ((3*_a3*(4*minus_eps13e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) - (_a1*(4*minus_eps13e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr)))/3 - (2*B*(((_a1*(4*minus_eps13e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) - (3*_a3*(4*minus_eps13e)*pow(minus_eps_tr,3.0))/(2*pow(minus_epse_I2,2.5)))*(minus_eps22e) + ((3*_a3*(4*minus_eps13e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) - (_a1*(4*minus_eps13e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr)))/3 + (B*(((_a1*(4*minus_eps13e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) - (3*_a3*(4*minus_eps13e)*pow(minus_eps_tr,3.0))/(2*pow(minus_epse_I2,2.5)))*(minus_eps33e) + ((3*_a3*(4*minus_eps13e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) - (_a1*(4*minus_eps13e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr)))/3 - ((B - 1)*((gamma_damaged*(4*minus_eps13e))/(2*pow(minus_epse_I2,0.5)) - (gamma_damaged*(4*minus_eps13e)*(minus_eps11e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5))))/3 + (2*(B - 1)*((gamma_damaged*(4*minus_eps13e))/(2*pow(minus_epse_I2,0.5)) - (gamma_damaged*(4*minus_eps13e)*(minus_eps22e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5))))/3 - ((B - 1)*((gamma_damaged*(4*minus_eps13e))/(2*pow(minus_epse_I2,0.5)) - (gamma_damaged*(4*minus_eps13e)*(minus_eps33e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5))))/3);
          jacobian(1,5) = -pow(B,_m1)*_C_g*_dt*_m2*pow((2*((2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5))*(minus_eps22e) + (lambda + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr))*(B - 1))/3 - (((2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5))*(minus_eps11e) + (lambda + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr))*(B - 1))/3 - (((2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5))*(minus_eps33e) + (lambda + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr))*(B - 1))/3 + (B*((2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5))*(minus_eps11e) - ((3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) - 2*_a2 + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr)))/3 - (2*B*((2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5))*(minus_eps22e) - ((3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) - 2*_a2 + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr)))/3 + (B*((2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5))*(minus_eps33e) - ((3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) - 2*_a2 + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr)))/3,_m2 - 1)*((B*(((_a1*(4*minus_eps23e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) - (3*_a3*(4*minus_eps23e)*pow(minus_eps_tr,3.0))/(2*pow(minus_epse_I2,2.5)))*(minus_eps11e) + ((3*_a3*(4*minus_eps23e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) - (_a1*(4*minus_eps23e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr)))/3 - (2*B*(((_a1*(4*minus_eps23e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) - (3*_a3*(4*minus_eps23e)*pow(minus_eps_tr,3.0))/(2*pow(minus_epse_I2,2.5)))*(minus_eps22e) + ((3*_a3*(4*minus_eps23e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) - (_a1*(4*minus_eps23e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr)))/3 + (B*(((_a1*(4*minus_eps23e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) - (3*_a3*(4*minus_eps23e)*pow(minus_eps_tr,3.0))/(2*pow(minus_epse_I2,2.5)))*(minus_eps33e) + ((3*_a3*(4*minus_eps23e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) - (_a1*(4*minus_eps23e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr)))/3 - ((B - 1)*((gamma_damaged*(4*minus_eps23e))/(2*pow(minus_epse_I2,0.5)) - (gamma_damaged*(4*minus_eps23e)*(minus_eps11e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5))))/3 + (2*(B - 1)*((gamma_damaged*(4*minus_eps23e))/(2*pow(minus_epse_I2,0.5)) - (gamma_damaged*(4*minus_eps23e)*(minus_eps22e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5))))/3 - ((B - 1)*((gamma_damaged*(4*minus_eps23e))/(2*pow(minus_epse_I2,0.5)) - (gamma_damaged*(4*minus_eps23e)*(minus_eps33e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5))))/3);
          jacobian(2,0) = pow(B,_m1)*_C_g*_dt*_m2*pow((2*((2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5))*(minus_eps33e) + (lambda + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr))*(B - 1))/3 - (((2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5))*(minus_eps22e) + (lambda + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr))*(B - 1))/3 - (((2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5))*(minus_eps11e) + (lambda + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr))*(B - 1))/3 + (B*((2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5))*(minus_eps11e) - ((3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) - 2*_a2 + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr)))/3 + (B*((2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5))*(minus_eps22e) - ((3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) - 2*_a2 + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr)))/3 - (2*B*((2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5))*(minus_eps33e) - ((3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) - 2*_a2 + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr)))/3,_m2 - 1)*((B*((minus_eps11e)*(_a1/pow(minus_epse_I2,0.5) - (3*_a3*pow(minus_eps_tr,2.0))/pow(minus_epse_I2,1.5) - (_a1*(2*minus_eps11e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) + (3*_a3*(2*minus_eps11e)*pow(minus_eps_tr,3.0))/(2*pow(minus_epse_I2,2.5))) - 2*_a2 - 2*_a0 + ((3*_a3)/pow(minus_epse_I2,0.5) - (_a1*pow(minus_epse_I2,0.5))/pow(minus_eps_tr,2.0) - (3*_a3*(2*minus_eps11e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) + (_a1*(2*minus_eps11e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr) + (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr) - (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5)))/3 + ((B - 1)*(lambda - ((gamma_damaged*pow(minus_epse_I2,0.5))/pow(minus_eps_tr,2.0) - (gamma_damaged*(2*minus_eps11e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr) + (gamma_damaged/pow(minus_epse_I2,0.5) - (gamma_damaged*(2*minus_eps11e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)))*(minus_eps22e) + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr)))/3 - (2*(B - 1)*(lambda - ((gamma_damaged*pow(minus_epse_I2,0.5))/pow(minus_eps_tr,2.0) - (gamma_damaged*(2*minus_eps11e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr) + (gamma_damaged/pow(minus_epse_I2,0.5) - (gamma_damaged*(2*minus_eps11e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)))*(minus_eps33e) + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr)))/3 + ((B - 1)*(lambda + 2*shear_modulus - ((gamma_damaged*pow(minus_epse_I2,0.5))/pow(minus_eps_tr,2.0) - (gamma_damaged*(2*minus_eps11e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr) + (gamma_damaged/pow(minus_epse_I2,0.5) - (gamma_damaged*(2*minus_eps11e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)))*(minus_eps11e) + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr)))/3 + (B*((minus_eps22e)*(_a1/pow(minus_epse_I2,0.5) - (3*_a3*pow(minus_eps_tr,2.0))/pow(minus_epse_I2,1.5) - (_a1*(2*minus_eps11e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) + (3*_a3*(2*minus_eps11e)*pow(minus_eps_tr,3.0))/(2*pow(minus_epse_I2,2.5))) - 2*_a2 + ((3*_a3)/pow(minus_epse_I2,0.5) - (_a1*pow(minus_epse_I2,0.5))/pow(minus_eps_tr,2.0) - (3*_a3*(2*minus_eps11e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) + (_a1*(2*minus_eps11e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr) + (3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr)))/3 - (2*B*((minus_eps33e)*(_a1/pow(minus_epse_I2,0.5) - (3*_a3*pow(minus_eps_tr,2.0))/pow(minus_epse_I2,1.5) - (_a1*(2*minus_eps11e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) + (3*_a3*(2*minus_eps11e)*pow(minus_eps_tr,3.0))/(2*pow(minus_epse_I2,2.5))) - 2*_a2 + ((3*_a3)/pow(minus_epse_I2,0.5) - (_a1*pow(minus_epse_I2,0.5))/pow(minus_eps_tr,2.0) - (3*_a3*(2*minus_eps11e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) + (_a1*(2*minus_eps11e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr) + (3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr)))/3);
          jacobian(2,1) = pow(B,_m1)*_C_g*_dt*_m2*pow((2*((2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5))*(minus_eps33e) + (lambda + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr))*(B - 1))/3 - (((2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5))*(minus_eps22e) + (lambda + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr))*(B - 1))/3 - (((2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5))*(minus_eps11e) + (lambda + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr))*(B - 1))/3 + (B*((2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5))*(minus_eps11e) - ((3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) - 2*_a2 + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr)))/3 + (B*((2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5))*(minus_eps22e) - ((3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) - 2*_a2 + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr)))/3 - (2*B*((2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5))*(minus_eps33e) - ((3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) - 2*_a2 + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr)))/3,_m2 - 1)*((B*((minus_eps22e)*(_a1/pow(minus_epse_I2,0.5) - (3*_a3*pow(minus_eps_tr,2.0))/pow(minus_epse_I2,1.5) - (_a1*(2*minus_eps22e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) + (3*_a3*(2*minus_eps22e)*pow(minus_eps_tr,3.0))/(2*pow(minus_epse_I2,2.5))) - 2*_a2 - 2*_a0 + ((3*_a3)/pow(minus_epse_I2,0.5) - (_a1*pow(minus_epse_I2,0.5))/pow(minus_eps_tr,2.0) - (3*_a3*(2*minus_eps22e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) + (_a1*(2*minus_eps22e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr) + (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr) - (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5)))/3 + ((B - 1)*(lambda - ((gamma_damaged*pow(minus_epse_I2,0.5))/pow(minus_eps_tr,2.0) - (gamma_damaged*(2*minus_eps22e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr) + (gamma_damaged/pow(minus_epse_I2,0.5) - (gamma_damaged*(2*minus_eps22e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)))*(minus_eps11e) + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr)))/3 - (2*(B - 1)*(lambda - ((gamma_damaged*pow(minus_epse_I2,0.5))/pow(minus_eps_tr,2.0) - (gamma_damaged*(2*minus_eps22e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr) + (gamma_damaged/pow(minus_epse_I2,0.5) - (gamma_damaged*(2*minus_eps22e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)))*(minus_eps33e) + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr)))/3 + ((B - 1)*(lambda + 2*shear_modulus - ((gamma_damaged*pow(minus_epse_I2,0.5))/pow(minus_eps_tr,2.0) - (gamma_damaged*(2*minus_eps22e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr) + (gamma_damaged/pow(minus_epse_I2,0.5) - (gamma_damaged*(2*minus_eps22e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)))*(minus_eps22e) + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr)))/3 + (B*((minus_eps11e)*(_a1/pow(minus_epse_I2,0.5) - (3*_a3*pow(minus_eps_tr,2.0))/pow(minus_epse_I2,1.5) - (_a1*(2*minus_eps22e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) + (3*_a3*(2*minus_eps22e)*pow(minus_eps_tr,3.0))/(2*pow(minus_epse_I2,2.5))) - 2*_a2 + ((3*_a3)/pow(minus_epse_I2,0.5) - (_a1*pow(minus_epse_I2,0.5))/pow(minus_eps_tr,2.0) - (3*_a3*(2*minus_eps22e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) + (_a1*(2*minus_eps22e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr) + (3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr)))/3 - (2*B*((minus_eps33e)*(_a1/pow(minus_epse_I2,0.5) - (3*_a3*pow(minus_eps_tr,2.0))/pow(minus_epse_I2,1.5) - (_a1*(2*minus_eps22e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) + (3*_a3*(2*minus_eps22e)*pow(minus_eps_tr,3.0))/(2*pow(minus_epse_I2,2.5))) - 2*_a2 + ((3*_a3)/pow(minus_epse_I2,0.5) - (_a1*pow(minus_epse_I2,0.5))/pow(minus_eps_tr,2.0) - (3*_a3*(2*minus_eps22e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) + (_a1*(2*minus_eps22e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr) + (3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr)))/3);
          jacobian(2,2) = pow(B,_m1)*_C_g*_dt*_m2*pow((2*((2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5))*(minus_eps33e) + (lambda + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr))*(B - 1))/3 - (((2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5))*(minus_eps22e) + (lambda + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr))*(B - 1))/3 - (((2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5))*(minus_eps11e) + (lambda + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr))*(B - 1))/3 + (B*((2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5))*(minus_eps11e) - ((3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) - 2*_a2 + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr)))/3 + (B*((2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5))*(minus_eps22e) - ((3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) - 2*_a2 + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr)))/3 - (2*B*((2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5))*(minus_eps33e) - ((3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) - 2*_a2 + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr)))/3,_m2 - 1)*(((B - 1)*(lambda - ((gamma_damaged*pow(minus_epse_I2,0.5))/pow(minus_eps_tr,2.0) - (gamma_damaged*(2*minus_eps33e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr) + (gamma_damaged/pow(minus_epse_I2,0.5) - (gamma_damaged*(2*minus_eps33e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)))*(minus_eps11e) + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr)))/3 - (2*B*((minus_eps33e)*(_a1/pow(minus_epse_I2,0.5) - (3*_a3*pow(minus_eps_tr,2.0))/pow(minus_epse_I2,1.5) - (_a1*(2*minus_eps33e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) + (3*_a3*(2*minus_eps33e)*pow(minus_eps_tr,3.0))/(2*pow(minus_epse_I2,2.5))) - 2*_a2 - 2*_a0 + ((3*_a3)/pow(minus_epse_I2,0.5) - (_a1*pow(minus_epse_I2,0.5))/pow(minus_eps_tr,2.0) - (3*_a3*(2*minus_eps33e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) + (_a1*(2*minus_eps33e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr) + (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr) - (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5)))/3 + ((B - 1)*(lambda - ((gamma_damaged*pow(minus_epse_I2,0.5))/pow(minus_eps_tr,2.0) - (gamma_damaged*(2*minus_eps33e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr) + (gamma_damaged/pow(minus_epse_I2,0.5) - (gamma_damaged*(2*minus_eps33e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)))*(minus_eps22e) + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr)))/3 - (2*(B - 1)*(lambda + 2*shear_modulus - ((gamma_damaged*pow(minus_epse_I2,0.5))/pow(minus_eps_tr,2.0) - (gamma_damaged*(2*minus_eps33e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr) + (gamma_damaged/pow(minus_epse_I2,0.5) - (gamma_damaged*(2*minus_eps33e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)))*(minus_eps33e) + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr)))/3 + (B*((minus_eps11e)*(_a1/pow(minus_epse_I2,0.5) - (3*_a3*pow(minus_eps_tr,2.0))/pow(minus_epse_I2,1.5) - (_a1*(2*minus_eps33e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) + (3*_a3*(2*minus_eps33e)*pow(minus_eps_tr,3.0))/(2*pow(minus_epse_I2,2.5))) - 2*_a2 + ((3*_a3)/pow(minus_epse_I2,0.5) - (_a1*pow(minus_epse_I2,0.5))/pow(minus_eps_tr,2.0) - (3*_a3*(2*minus_eps33e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) + (_a1*(2*minus_eps33e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr) + (3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr)))/3 + (B*((minus_eps22e)*(_a1/pow(minus_epse_I2,0.5) - (3*_a3*pow(minus_eps_tr,2.0))/pow(minus_epse_I2,1.5) - (_a1*(2*minus_eps33e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) + (3*_a3*(2*minus_eps33e)*pow(minus_eps_tr,3.0))/(2*pow(minus_epse_I2,2.5))) - 2*_a2 + ((3*_a3)/pow(minus_epse_I2,0.5) - (_a1*pow(minus_epse_I2,0.5))/pow(minus_eps_tr,2.0) - (3*_a3*(2*minus_eps33e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) + (_a1*(2*minus_eps33e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr) + (3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr)))/3) + 1;
          jacobian(2,3) = -pow(B,_m1)*_C_g*_dt*_m2*pow((2*((2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5))*(minus_eps33e) + (lambda + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr))*(B - 1))/3 - (((2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5))*(minus_eps22e) + (lambda + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr))*(B - 1))/3 - (((2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5))*(minus_eps11e) + (lambda + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr))*(B - 1))/3 + (B*((2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5))*(minus_eps11e) - ((3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) - 2*_a2 + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr)))/3 + (B*((2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5))*(minus_eps22e) - ((3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) - 2*_a2 + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr)))/3 - (2*B*((2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5))*(minus_eps33e) - ((3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) - 2*_a2 + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr)))/3,_m2 - 1)*((B*(((_a1*(4*minus_eps12e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) - (3*_a3*(4*minus_eps12e)*pow(minus_eps_tr,3.0))/(2*pow(minus_epse_I2,2.5)))*(minus_eps11e) + ((3*_a3*(4*minus_eps12e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) - (_a1*(4*minus_eps12e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr)))/3 + (B*(((_a1*(4*minus_eps12e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) - (3*_a3*(4*minus_eps12e)*pow(minus_eps_tr,3.0))/(2*pow(minus_epse_I2,2.5)))*(minus_eps22e) + ((3*_a3*(4*minus_eps12e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) - (_a1*(4*minus_eps12e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr)))/3 - (2*B*(((_a1*(4*minus_eps12e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) - (3*_a3*(4*minus_eps12e)*pow(minus_eps_tr,3.0))/(2*pow(minus_epse_I2,2.5)))*(minus_eps33e) + ((3*_a3*(4*minus_eps12e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) - (_a1*(4*minus_eps12e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr)))/3 - ((B - 1)*((gamma_damaged*(4*minus_eps12e))/(2*pow(minus_epse_I2,0.5)) - (gamma_damaged*(4*minus_eps12e)*(minus_eps11e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5))))/3 - ((B - 1)*((gamma_damaged*(4*minus_eps12e))/(2*pow(minus_epse_I2,0.5)) - (gamma_damaged*(4*minus_eps12e)*(minus_eps22e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5))))/3 + (2*(B - 1)*((gamma_damaged*(4*minus_eps12e))/(2*pow(minus_epse_I2,0.5)) - (gamma_damaged*(4*minus_eps12e)*(minus_eps33e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5))))/3);
          jacobian(2,4) = -pow(B,_m1)*_C_g*_dt*_m2*pow((2*((2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5))*(minus_eps33e) + (lambda + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr))*(B - 1))/3 - (((2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5))*(minus_eps22e) + (lambda + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr))*(B - 1))/3 - (((2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5))*(minus_eps11e) + (lambda + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr))*(B - 1))/3 + (B*((2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5))*(minus_eps11e) - ((3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) - 2*_a2 + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr)))/3 + (B*((2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5))*(minus_eps22e) - ((3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) - 2*_a2 + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr)))/3 - (2*B*((2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5))*(minus_eps33e) - ((3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) - 2*_a2 + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr)))/3,_m2 - 1)*((B*(((_a1*(4*minus_eps13e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) - (3*_a3*(4*minus_eps13e)*pow(minus_eps_tr,3.0))/(2*pow(minus_epse_I2,2.5)))*(minus_eps11e) + ((3*_a3*(4*minus_eps13e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) - (_a1*(4*minus_eps13e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr)))/3 + (B*(((_a1*(4*minus_eps13e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) - (3*_a3*(4*minus_eps13e)*pow(minus_eps_tr,3.0))/(2*pow(minus_epse_I2,2.5)))*(minus_eps22e) + ((3*_a3*(4*minus_eps13e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) - (_a1*(4*minus_eps13e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr)))/3 - (2*B*(((_a1*(4*minus_eps13e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) - (3*_a3*(4*minus_eps13e)*pow(minus_eps_tr,3.0))/(2*pow(minus_epse_I2,2.5)))*(minus_eps33e) + ((3*_a3*(4*minus_eps13e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) - (_a1*(4*minus_eps13e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr)))/3 - ((B - 1)*((gamma_damaged*(4*minus_eps13e))/(2*pow(minus_epse_I2,0.5)) - (gamma_damaged*(4*minus_eps13e)*(minus_eps11e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5))))/3 - ((B - 1)*((gamma_damaged*(4*minus_eps13e))/(2*pow(minus_epse_I2,0.5)) - (gamma_damaged*(4*minus_eps13e)*(minus_eps22e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5))))/3 + (2*(B - 1)*((gamma_damaged*(4*minus_eps13e))/(2*pow(minus_epse_I2,0.5)) - (gamma_damaged*(4*minus_eps13e)*(minus_eps33e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5))))/3);
          jacobian(2,5) = -pow(B,_m1)*_C_g*_dt*_m2*pow((2*((2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5))*(minus_eps33e) + (lambda + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr))*(B - 1))/3 - (((2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5))*(minus_eps22e) + (lambda + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr))*(B - 1))/3 - (((2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5))*(minus_eps11e) + (lambda + (gamma_damaged*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr))*(B - 1))/3 + (B*((2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5))*(minus_eps11e) - ((3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) - 2*_a2 + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr)))/3 + (B*((2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5))*(minus_eps22e) - ((3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) - 2*_a2 + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr)))/3 - (2*B*((2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5))*(minus_eps33e) - ((3*_a3*(minus_eps_tr))/pow(minus_epse_I2,0.5) - 2*_a2 + (_a1*pow(minus_epse_I2,0.5))/(minus_eps_tr))*(minus_eps_tr)))/3,_m2 - 1)*((B*(((_a1*(4*minus_eps23e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) - (3*_a3*(4*minus_eps23e)*pow(minus_eps_tr,3.0))/(2*pow(minus_epse_I2,2.5)))*(minus_eps11e) + ((3*_a3*(4*minus_eps23e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) - (_a1*(4*minus_eps23e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr)))/3 + (B*(((_a1*(4*minus_eps23e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) - (3*_a3*(4*minus_eps23e)*pow(minus_eps_tr,3.0))/(2*pow(minus_epse_I2,2.5)))*(minus_eps22e) + ((3*_a3*(4*minus_eps23e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) - (_a1*(4*minus_eps23e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr)))/3 - (2*B*(((_a1*(4*minus_eps23e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) - (3*_a3*(4*minus_eps23e)*pow(minus_eps_tr,3.0))/(2*pow(minus_epse_I2,2.5)))*(minus_eps33e) + ((3*_a3*(4*minus_eps23e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) - (_a1*(4*minus_eps23e))/(2*pow(minus_epse_I2,0.5)*(minus_eps_tr)))*(minus_eps_tr)))/3 - ((B - 1)*((gamma_damaged*(4*minus_eps23e))/(2*pow(minus_epse_I2,0.5)) - (gamma_damaged*(4*minus_eps23e)*(minus_eps11e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5))))/3 - ((B - 1)*((gamma_damaged*(4*minus_eps23e))/(2*pow(minus_epse_I2,0.5)) - (gamma_damaged*(4*minus_eps23e)*(minus_eps22e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5))))/3 + (2*(B - 1)*((gamma_damaged*(4*minus_eps23e))/(2*pow(minus_epse_I2,0.5)) - (gamma_damaged*(4*minus_eps23e)*(minus_eps33e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5))))/3);
          jacobian(3,0) = -pow(B,_m1)*_C_g*_dt*_m2*pow((B - 1)*(2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5))*(minus_eps12e) - B*(2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5))*(minus_eps12e),_m2 - 1)*((gamma_damaged/pow(minus_epse_I2,0.5) - (gamma_damaged*(2*minus_eps11e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)))*(B - 1)*(minus_eps12e) + B*(minus_eps12e)*(_a1/pow(minus_epse_I2,0.5) - (3*_a3*pow(minus_eps_tr,2.0))/pow(minus_epse_I2,1.5) - (_a1*(2*minus_eps11e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) + (3*_a3*(2*minus_eps11e)*pow(minus_eps_tr,3.0))/(2*pow(minus_epse_I2,2.5))));
          jacobian(3,1) = -pow(B,_m1)*_C_g*_dt*_m2*pow((B - 1)*(2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5))*(minus_eps12e) - B*(2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5))*(minus_eps12e),_m2 - 1)*((gamma_damaged/pow(minus_epse_I2,0.5) - (gamma_damaged*(2*minus_eps22e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)))*(B - 1)*(minus_eps12e) + B*(minus_eps12e)*(_a1/pow(minus_epse_I2,0.5) - (3*_a3*pow(minus_eps_tr,2.0))/pow(minus_epse_I2,1.5) - (_a1*(2*minus_eps22e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) + (3*_a3*(2*minus_eps22e)*pow(minus_eps_tr,3.0))/(2*pow(minus_epse_I2,2.5))));
          jacobian(3,2) = -pow(B,_m1)*_C_g*_dt*_m2*pow((B - 1)*(2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5))*(minus_eps12e) - B*(2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5))*(minus_eps12e),_m2 - 1)*((gamma_damaged/pow(minus_epse_I2,0.5) - (gamma_damaged*(2*minus_eps33e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)))*(B - 1)*(minus_eps12e) + B*(minus_eps12e)*(_a1/pow(minus_epse_I2,0.5) - (3*_a3*pow(minus_eps_tr,2.0))/pow(minus_epse_I2,1.5) - (_a1*(2*minus_eps33e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) + (3*_a3*(2*minus_eps33e)*pow(minus_eps_tr,3.0))/(2*pow(minus_epse_I2,2.5))));
          jacobian(3,3) = pow(B,_m1)*_C_g*_dt*_m2*pow((B - 1)*(2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5))*(minus_eps12e) - B*(2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5))*(minus_eps12e),_m2 - 1)*(B*(2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5)) - (B - 1)*(2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5)) + B*((_a1*(4*minus_eps12e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) - (3*_a3*(4*minus_eps12e)*pow(minus_eps_tr,3.0))/(2*pow(minus_epse_I2,2.5)))*(minus_eps12e) + (gamma_damaged*(B - 1)*(4*minus_eps12e)*(minus_eps12e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5))) + 1;
          jacobian(3,4) = pow(B,_m1)*_C_g*_dt*_m2*pow((B - 1)*(2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5))*(minus_eps12e) - B*(2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5))*(minus_eps12e),_m2 - 1)*(B*((_a1*(4*minus_eps13e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) - (3*_a3*(4*minus_eps13e)*pow(minus_eps_tr,3.0))/(2*pow(minus_epse_I2,2.5)))*(minus_eps12e) + (gamma_damaged*(B - 1)*(4*minus_eps13e)*(minus_eps12e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)));
          jacobian(3,5) = pow(B,_m1)*_C_g*_dt*_m2*pow((B - 1)*(2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5))*(minus_eps12e) - B*(2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5))*(minus_eps12e),_m2 - 1)*(B*((_a1*(4*minus_eps23e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) - (3*_a3*(4*minus_eps23e)*pow(minus_eps_tr,3.0))/(2*pow(minus_epse_I2,2.5)))*(minus_eps12e) + (gamma_damaged*(B - 1)*(4*minus_eps23e)*(minus_eps12e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5))); 
          jacobian(4,0) = -pow(B,_m1)*_C_g*_dt*_m2*pow((B - 1)*(2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5))*(minus_eps13e) - B*(2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5))*(minus_eps13e),_m2 - 1)*((gamma_damaged/pow(minus_epse_I2,0.5) - (gamma_damaged*(2*minus_eps11e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)))*(B - 1)*(minus_eps13e) + B*(minus_eps13e)*(_a1/pow(minus_epse_I2,0.5) - (3*_a3*pow(minus_eps_tr,2.0))/pow(minus_epse_I2,1.5) - (_a1*(2*minus_eps11e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) + (3*_a3*(2*minus_eps11e)*pow(minus_eps_tr,3.0))/(2*pow(minus_epse_I2,2.5))));
          jacobian(4,1) = -pow(B,_m1)*_C_g*_dt*_m2*pow((B - 1)*(2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5))*(minus_eps13e) - B*(2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5))*(minus_eps13e),_m2 - 1)*((gamma_damaged/pow(minus_epse_I2,0.5) - (gamma_damaged*(2*minus_eps22e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)))*(B - 1)*(minus_eps13e) + B*(minus_eps13e)*(_a1/pow(minus_epse_I2,0.5) - (3*_a3*pow(minus_eps_tr,2.0))/pow(minus_epse_I2,1.5) - (_a1*(2*minus_eps22e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) + (3*_a3*(2*minus_eps22e)*pow(minus_eps_tr,3.0))/(2*pow(minus_epse_I2,2.5))));
          jacobian(4,2) = -pow(B,_m1)*_C_g*_dt*_m2*pow((B - 1)*(2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5))*(minus_eps13e) - B*(2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5))*(minus_eps13e),_m2 - 1)*((gamma_damaged/pow(minus_epse_I2,0.5) - (gamma_damaged*(2*minus_eps33e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)))*(B - 1)*(minus_eps13e) + B*(minus_eps13e)*(_a1/pow(minus_epse_I2,0.5) - (3*_a3*pow(minus_eps_tr,2.0))/pow(minus_epse_I2,1.5) - (_a1*(2*minus_eps33e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) + (3*_a3*(2*minus_eps33e)*pow(minus_eps_tr,3.0))/(2*pow(minus_epse_I2,2.5))));
          jacobian(4,3) = pow(B,_m1)*_C_g*_dt*_m2*pow((B - 1)*(2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5))*(minus_eps13e) - B*(2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5))*(minus_eps13e),_m2 - 1)*(B*((_a1*(4*minus_eps12e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) - (3*_a3*(4*minus_eps12e)*pow(minus_eps_tr,3.0))/(2*pow(minus_epse_I2,2.5)))*(minus_eps13e) + (gamma_damaged*(B - 1)*(4*minus_eps12e)*(minus_eps13e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)));
          jacobian(4,4) = pow(B,_m1)*_C_g*_dt*_m2*pow((B - 1)*(2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5))*(minus_eps13e) - B*(2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5))*(minus_eps13e),_m2 - 1)*(B*(2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5)) - (B - 1)*(2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5)) + B*((_a1*(4*minus_eps13e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) - (3*_a3*(4*minus_eps13e)*pow(minus_eps_tr,3.0))/(2*pow(minus_epse_I2,2.5)))*(minus_eps13e) + (gamma_damaged*(B - 1)*(4*minus_eps13e)*(minus_eps13e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5))) + 1;
          jacobian(4,5) = pow(B,_m1)*_C_g*_dt*_m2*pow((B - 1)*(2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5))*(minus_eps13e) - B*(2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5))*(minus_eps13e),_m2 - 1)*(B*((_a1*(4*minus_eps23e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) - (3*_a3*(4*minus_eps23e)*pow(minus_eps_tr,3.0))/(2*pow(minus_epse_I2,2.5)))*(minus_eps13e) + (gamma_damaged*(B - 1)*(4*minus_eps23e)*(minus_eps13e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5))); 
          jacobian(5,0) = -pow(B,_m1)*_C_g*_dt*_m2*pow((B - 1)*(2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5))*(minus_eps23e) - B*(2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5))*(minus_eps23e),_m2 - 1)*((gamma_damaged/pow(minus_epse_I2,0.5) - (gamma_damaged*(2*minus_eps11e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)))*(B - 1)*(minus_eps23e) + B*(minus_eps23e)*(_a1/pow(minus_epse_I2,0.5) - (3*_a3*pow(minus_eps_tr,2.0))/pow(minus_epse_I2,1.5) - (_a1*(2*minus_eps11e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) + (3*_a3*(2*minus_eps11e)*pow(minus_eps_tr,3.0))/(2*pow(minus_epse_I2,2.5))));
          jacobian(5,1) = -pow(B,_m1)*_C_g*_dt*_m2*pow((B - 1)*(2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5))*(minus_eps23e) - B*(2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5))*(minus_eps23e),_m2 - 1)*((gamma_damaged/pow(minus_epse_I2,0.5) - (gamma_damaged*(2*minus_eps22e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)))*(B - 1)*(minus_eps23e) + B*(minus_eps23e)*(_a1/pow(minus_epse_I2,0.5) - (3*_a3*pow(minus_eps_tr,2.0))/pow(minus_epse_I2,1.5) - (_a1*(2*minus_eps22e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) + (3*_a3*(2*minus_eps22e)*pow(minus_eps_tr,3.0))/(2*pow(minus_epse_I2,2.5))));
          jacobian(5,2) = -pow(B,_m1)*_C_g*_dt*_m2*pow((B - 1)*(2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5))*(minus_eps23e) - B*(2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5))*(minus_eps23e),_m2 - 1)*((gamma_damaged/pow(minus_epse_I2,0.5) - (gamma_damaged*(2*minus_eps33e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)))*(B - 1)*(minus_eps23e) + B*(minus_eps23e)*(_a1/pow(minus_epse_I2,0.5) - (3*_a3*pow(minus_eps_tr,2.0))/pow(minus_epse_I2,1.5) - (_a1*(2*minus_eps33e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) + (3*_a3*(2*minus_eps33e)*pow(minus_eps_tr,3.0))/(2*pow(minus_epse_I2,2.5))));
          jacobian(5,3) = pow(B,_m1)*_C_g*_dt*_m2*pow((B - 1)*(2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5))*(minus_eps23e) - B*(2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5))*(minus_eps23e),_m2 - 1)*(B*((_a1*(4*minus_eps12e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) - (3*_a3*(4*minus_eps12e)*pow(minus_eps_tr,3.0))/(2*pow(minus_epse_I2,2.5)))*(minus_eps23e) + (gamma_damaged*(B - 1)*(4*minus_eps12e)*(minus_eps23e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)));
          jacobian(5,4) = pow(B,_m1)*_C_g*_dt*_m2*pow((B - 1)*(2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5))*(minus_eps23e) - B*(2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5))*(minus_eps23e),_m2 - 1)*(B*((_a1*(4*minus_eps13e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) - (3*_a3*(4*minus_eps13e)*pow(minus_eps_tr,3.0))/(2*pow(minus_epse_I2,2.5)))*(minus_eps23e) + (gamma_damaged*(B - 1)*(4*minus_eps13e)*(minus_eps23e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)));
          jacobian(5,5) = pow(B,_m1)*_C_g*_dt*_m2*pow((B - 1)*(2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5))*(minus_eps23e) - B*(2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5))*(minus_eps23e),_m2 - 1)*(B*(2*_a0 - (_a1*(minus_eps_tr))/pow(minus_epse_I2,0.5) + (_a3*pow(minus_eps_tr,3.0))/pow(minus_epse_I2,1.5)) - (B - 1)*(2*shear_modulus + (gamma_damaged*(minus_eps_tr))/pow(minus_epse_I2,0.5)) + B*((_a1*(4*minus_eps23e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5)) - (3*_a3*(4*minus_eps23e)*pow(minus_eps_tr,3.0))/(2*pow(minus_epse_I2,2.5)))*(minus_eps23e) + (gamma_damaged*(B - 1)*(4*minus_eps23e)*(minus_eps23e)*(minus_eps_tr))/(2*pow(minus_epse_I2,1.5))) + 1; 
      };   

      NestedSolve solver;
      NestedSolve::Value<> solution(6);
      solution << 0, 0, 0, 0, 0, 0;
      solver.setAbsoluteTolerance(1e-10);
      solver.setRelativeTolerance(1e-8);
      solver.nonlinear(solution, compute);

      //take solution
      Real eps11p_inc = solution(0);
      Real eps22p_inc = solution(1);
      Real eps33p_inc = solution(2);
      Real eps12p_inc = solution(3);
      Real eps13p_inc = solution(4);
      Real eps23p_inc = solution(5);

      //save data
      //viscoelastic strain
      RankTwoTensor eps_p_out;
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
      RankTwoTensor eps_e_out;
      Real eps11e_out = eps11t - ( eps11p_pre + eps11p_inc );
      Real eps22e_out = eps22t - ( eps22p_pre + eps22p_inc );
      Real eps33e_out = eps33t - ( eps33p_pre + eps33p_inc );
      Real eps12e_out = eps12t - ( eps12p_pre + eps12p_inc );
      Real eps13e_out = eps13t - ( eps13p_pre + eps13p_inc );
      Real eps23e_out = eps23t - ( eps23p_pre + eps23p_inc );
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
      RankTwoTensor eps_total_out;
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
      
      Real I1_out = eps11e_out + eps22e_out + eps33e_out;
      Real I2_out = eps11e_out * eps11e_out + eps22e_out * eps22e_out + eps33e_out * eps33e_out + 2 * eps12e_out * eps12e_out + 2 * eps13e_out * eps13e_out + 2 * eps23e_out * eps23e_out;
      Real xi_out = I1_out / sqrt(I2_out);

      //save
      _I1[_qp] = I1_out;
      _I2[_qp] = I2_out;
      _xi[_qp] = xi_out; 

      //compute stress
      //sts_total, stress are updated
      //feed total stress
      RankTwoTensor stress_total_out;
      stress_total_out(0,0) = computeStressComps(1, 1, xi_out, I1_out, B_out, lambda_out, gamma_damaged_out, shear_modulus_out, eps11e_out, eps22e_out, eps12e_out, eps33e_out, eps13e_out, eps23e_out, alpha_grad_x, alpha_grad_y, alpha_grad_z, D);
      stress_total_out(1,1) = computeStressComps(2, 2, xi_out, I1_out, B_out, lambda_out, gamma_damaged_out, shear_modulus_out, eps11e_out, eps22e_out, eps12e_out, eps33e_out, eps13e_out, eps23e_out, alpha_grad_x, alpha_grad_y, alpha_grad_z, D);
      stress_total_out(0,1) = computeStressComps(1, 2, xi_out, I1_out, B_out, lambda_out, gamma_damaged_out, shear_modulus_out, eps11e_out, eps22e_out, eps12e_out, eps33e_out, eps13e_out, eps23e_out, alpha_grad_x, alpha_grad_y, alpha_grad_z, D);
      stress_total_out(1,0) = computeStressComps(1, 2, xi_out, I1_out, B_out, lambda_out, gamma_damaged_out, shear_modulus_out, eps11e_out, eps22e_out, eps12e_out, eps33e_out, eps13e_out, eps23e_out, alpha_grad_x, alpha_grad_y, alpha_grad_z, D);
      stress_total_out(2,2) = computeStressComps(3, 3, xi_out, I1_out, B_out, lambda_out, gamma_damaged_out, shear_modulus_out, eps11e_out, eps22e_out, eps12e_out, eps33e_out, eps13e_out, eps23e_out, alpha_grad_x, alpha_grad_y, alpha_grad_z, D);
      stress_total_out(0,2) = computeStressComps(1, 3, xi_out, I1_out, B_out, lambda_out, gamma_damaged_out, shear_modulus_out, eps11e_out, eps22e_out, eps12e_out, eps33e_out, eps13e_out, eps23e_out, alpha_grad_x, alpha_grad_y, alpha_grad_z, D);
      stress_total_out(2,0) = computeStressComps(1, 3, xi_out, I1_out, B_out, lambda_out, gamma_damaged_out, shear_modulus_out, eps11e_out, eps22e_out, eps12e_out, eps33e_out, eps13e_out, eps23e_out, alpha_grad_x, alpha_grad_y, alpha_grad_z, D);
      stress_total_out(1,2) = computeStressComps(2, 3, xi_out, I1_out, B_out, lambda_out, gamma_damaged_out, shear_modulus_out, eps11e_out, eps22e_out, eps12e_out, eps33e_out, eps13e_out, eps23e_out, alpha_grad_x, alpha_grad_y, alpha_grad_z, D);
      stress_total_out(2,1) = computeStressComps(2, 3, xi_out, I1_out, B_out, lambda_out, gamma_damaged_out, shear_modulus_out, eps11e_out, eps22e_out, eps12e_out, eps33e_out, eps13e_out, eps23e_out, alpha_grad_x, alpha_grad_y, alpha_grad_z, D);

      _sts_total[_qp] = stress_total_out;

      //feed stress change (relative to initial condition) to system
      RankTwoTensor stress_out;
      RankTwoTensor stress_initial = _static_initial_stress_tensor[_qp];
      stress_out(0,0) = stress_total_out(0,0) - stress_initial(0,0);
      stress_out(1,1) = stress_total_out(1,1) - stress_initial(1,1);
      stress_out(2,2) = stress_total_out(2,2) - stress_initial(2,2);
      stress_out(0,1) = stress_total_out(0,1) - stress_initial(0,1);
      stress_out(1,0) = stress_total_out(1,0) - stress_initial(1,0);
      stress_out(0,2) = stress_total_out(0,2) - stress_initial(0,2);
      stress_out(2,0) = stress_total_out(2,0) - stress_initial(2,0);
      stress_out(1,2) = stress_total_out(1,2) - stress_initial(1,2);
      stress_out(2,1) = stress_total_out(2,1) - stress_initial(2,1);

      _stress[_qp] = stress_out; //this needs to feed in stress increment relative to initial value

      // stress = C * e
      //_stress[_qp] = _elasticity_tensor[_qp] * _mechanical_strain[_qp];

      // Assign value for elastic strain, which is equal to the mechanical strain
      _elastic_strain[_qp] = _mechanical_strain[_qp];

      // Compute dstress_dstrain (must use explicit solver)
      // _Jacobian_mult[_qp] = _elasticity_tensor[_qp];
    
    //}

  }
}

///Function: deltaij
Real 
ComputeDamageBreakageStress3D::deltaij(int i, int j)
{
  Real deltaij_out = 0.0;
  if ( i == j )
  {
    deltaij_out = 1.0;
  }
  else
  {
    deltaij_out = 0.0;
  }
  return deltaij_out;
}

/// Function: epsilonij - take component of elastic strain
Real 
ComputeDamageBreakageStress3D::epsilonij(int i, 
                                       int j,
                                       Real eps11e_in,
                                       Real eps22e_in,
                                       Real eps12e_in,
                                       Real eps33e_in,
                                       Real eps13e_in,
                                       Real eps23e_in)
{ 
  //take elastic strain
  Real eps_e_out = 0.0;
  if ( i == 1 && j == 1 )
  {
    eps_e_out = eps11e_in;
  }
  else if ( i == 2 && j == 2 )
  {
    eps_e_out = eps22e_in;
  }
  else if ( i == 3 && j == 3 )
  {
    eps_e_out = eps33e_in;
  }
  else if ( i == 1 && j == 2 )
  {
    eps_e_out = eps12e_in;
  }
  else if ( i == 1 && j == 3 )
  {
    eps_e_out = eps13e_in;
  }
  else if ( i == 2 && j == 3 )
  {
    eps_e_out = eps23e_in;
  }
  return eps_e_out;
}

/// Function: grad_alpha
//Note: no grad in 33 direction
Real 
ComputeDamageBreakageStress3D::grad_alpha(int i, 
                                          Real alpha_grad_x,
                                          Real alpha_grad_y,
                                          Real alpha_grad_z)
{ 
  //take alpha grad
  Real alpha_grad_out = 0.0;
  if ( i == 1 )
  {
    alpha_grad_out = alpha_grad_x;
  }
  else if ( i == 2 )
  {
    alpha_grad_out = alpha_grad_y;
  }
  else if ( i == 3 )
  {
    alpha_grad_out = alpha_grad_z;
  }
  return alpha_grad_out;
}

/// Function: compute stress components
Real 
ComputeDamageBreakageStress3D::computeStressComps(int i, 
                                                int j,
                                                Real xi_in,
                                                Real I1_in,
                                                Real B_in,
                                                Real lambda_in,
                                                Real gamma_damaged_in,
                                                Real shear_modulus_in,
                                                Real eps11e_in,
                                                Real eps22e_in,
                                                Real eps12e_in,
                                                Real eps33e_in,
                                                Real eps13e_in,
                                                Real eps23e_in,
                                                Real alpha_grad_x,
                                                Real alpha_grad_y,
                                                Real alpha_grad_z,
                                                Real D)
{
  //Retrieve parameters
  Real xi = xi_in;
  Real I1 = I1_in;
  Real B = B_in;
  Real lambda = lambda_in;
  Real gamma_damaged = gamma_damaged_in;
  Real shear_modulus = shear_modulus_in;

  //stress comps
  Real stresscomp_s_out;
  Real stresscomp_b_out;
  Real stresscomp;

  stresscomp_s_out = (lambda - gamma_damaged / xi) * I1 * deltaij(i,j) + ( 2 * shear_modulus - gamma_damaged * xi ) * epsilonij(i,j,eps11e_in,eps22e_in,eps12e_in,eps33e_in,eps13e_in,eps23e_in) - D * grad_alpha(i,alpha_grad_x,alpha_grad_y,alpha_grad_z) * grad_alpha(j,alpha_grad_x,alpha_grad_y,alpha_grad_z);
  stresscomp_b_out = (2 * _a2 + _a1 / xi + 3 * _a3 * xi) * I1 * deltaij(i,j) + ( 2 * _a0 + _a1 * xi - _a3 * pow(xi,3) ) * epsilonij(i,j,eps11e_in,eps22e_in,eps12e_in,eps33e_in,eps13e_in,eps23e_in);

  //
  stresscomp = ( 1 - B ) * stresscomp_s_out + B * stresscomp_b_out;

  return stresscomp;

}

///Function to compute initial strain based on initial stress 
void
ComputeDamageBreakageStress3D::setupInitial()
{
  ///For isotropic material with all components of stress subject to small strain we consider 
  ///tensile/compressive stress leads to only tensile/compressive strain, shear stress produce
  ///shear strain: 
  ///eps_ii = 1/E * ( sigma_ii - nu * ( sigma_jj + sigma_kk) )
  ///eps_ij = 1/G * sigma_ij

  /// lambda (first lame const)
  _lambda[_qp] = _lambda_o;
  /// mu (shear modulus)
  _shear_modulus[_qp] = _shear_modulus_o;
  /// gamma_damaged (damage modulus)
  _gamma_damaged[_qp] = 0.0;

  //allpha, B
  _alpha_damagedvar[_qp] = 0.0;
  _B[_qp] = 0.0;

  //Convert (lambda_o,shear_modulus_o) to (youngs_modulus_o,poisson_ratio_o)
  // Real youngs_modulus_o = _shear_modulus_o * ( 3 * _lambda_o + 2 * _shear_modulus_o ) / ( _lambda_o + _shear_modulus_o );
  // Real poisson_ratio_o = _lambda_o / ( 2 * ( _lambda_o + _shear_modulus_o ));

  //Convert (lambda_o,shear_modulus_o) to (shear_wave_speed_o,pressure_wave_speed_o)
  // Real density_o = _density_old[_qp];
  // Real shear_wave_speed_o = sqrt( ( _shear_modulus_o ) / ( density_o ) );
  // Real pressure_wave_speed_o = sqrt( ( _lambda_o + 2 * _shear_modulus_o ) / ( density_o ) );

  //save
  // _shear_wave_speed[_qp] = shear_wave_speed_o;
  // _pressure_wave_speed[_qp] = pressure_wave_speed_o;

  //Get stress components
  RankTwoTensor stress_initial = _static_initial_stress_tensor[_qp];

  //Get strain components
  RankTwoTensor strain_initial = _static_initial_strain_tensor[_qp];
  Real eps11_init = strain_initial(0,0);
  Real eps12_init = strain_initial(0,1);
  Real eps22_init = strain_initial(1,1);
  Real eps13_init = strain_initial(0,2);
  Real eps23_init = strain_initial(1,2);
  Real eps33_init = strain_initial(2,2);  

  //-------------------------------------------------------------------------------------------------------------------------------//
  //Initially we are at elastic region, eps33_init = eps_elastic = - eps_plastic = 0                                               //
  //So below we only compute 11 22 12 components, but notice the solving procedue should include out-of-plane strain components 33 //
  //when visco-elastic strain is invoked                                                                                           //
  //-------------------------------------------------------------------------------------------------------------------------------//

  //Compute xi, I1, I2
  Real I1_init = eps11_init + eps22_init + eps33_init;
  Real I2_init = eps11_init * eps11_init + eps22_init * eps22_init + eps33_init * eps33_init + 2 * eps12_init * eps12_init + 2 * eps13_init * eps13_init + 2 * eps23_init * eps23_init;
  Real xi_init = I1_init / sqrt( I2_init );

  //Compute eps
  //eps_p
  _eps_p[_qp](0,0) = 0.0; _eps_p[_qp](0,1) = 0.0; _eps_p[_qp](0,2) = 0.0;
  _eps_p[_qp](1,0) = 0.0; _eps_p[_qp](1,1) = 0.0; _eps_p[_qp](1,2) = 0.0;
  _eps_p[_qp](2,0) = 0.0; _eps_p[_qp](2,1) = 0.0; _eps_p[_qp](2,2) = 0.0;
  //eps_e
  _eps_e[_qp] = strain_initial;
  //eps_total
  _eps_total[_qp] = strain_initial;
  //sts_total
  _sts_total[_qp] = stress_initial;

  //I1
  _I1[_qp] = I1_init;
  //I2
  _I2[_qp] = I2_init;
  //xi
  _xi[_qp] = xi_init;
}