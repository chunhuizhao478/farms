//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DiffusedDamageBreakageMaterialMainApp.h"

/**
 *  Material used in damage-breakage large deformation formulation, consider full damage evolution equation with diffusion
 *  Created by Chunhui Zhao, Dec 24th, 2024
 */
registerMooseObject("farmsApp", DiffusedDamageBreakageMaterialMainApp);

InputParameters
DiffusedDamageBreakageMaterialMainApp::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Material used in three field poro dynamics simulations");
  //input parameters
  params.addRequiredParam<Real>(        "lambda_o", "initial lambda constant value");
  params.addRequiredParam<Real>( "shear_modulus_o", "initial shear modulus value");
  params.addRequiredParam<Real>(            "xi_0", "strain invariants ratio: onset of damage evolution");
  params.addRequiredParam<Real>(            "xi_d", "strain invariants ratio: onset of breakage healing");
  params.addRequiredParam<Real>(             "chi", "coefficient of energy ratio Fb/Fs = chi < 1");
  params.addRequiredParam<Real>(             "C_g", "compliance or fluidity of the fine grain granular material");
  params.addRequiredParam<Real>(              "m1", "coefficient of power law indexes");
  params.addRequiredParam<Real>(              "m2", "coefficient of power law indexes");  
  //input coupled variables from main app
  params.addRequiredCoupledVar("structural_stress_coefficient", "structral_stress_coefficient");
  params.addRequiredCoupledVar("alpha_damagedvar_aux", "second_elastic_strain_invariant");
  params.addRequiredCoupledVar("B_damagedvar_aux", "strain_invariant_ratio");
  //build L matrix
  params.addRequiredCoupledVar(           "vel_x", "velocity in x direction"); //to build L matrix
  params.addRequiredCoupledVar(           "vel_y", "velocity in y direction"); //to build L matrix
  params.addRequiredCoupledVar(           "vel_z", "velocity in z direction"); //to build L matrix  
  return params;
}

DiffusedDamageBreakageMaterialMainApp::DiffusedDamageBreakageMaterialMainApp(const InputParameters & parameters)
  : Material(parameters),
  //declare properties
  //--------------------------------------------------------------//
  _gamma_damaged_r(declareProperty<Real>("gamma_damaged_r")),
  _alpha_damagedvar(declareProperty<Real>("alpha_damagedvar")),
  _B_damagedvar(declareProperty<Real>("B_damagedvar")),
  _lambda(declareProperty<Real>("lambda_const")),
  _shear_modulus(declareProperty<Real>("shear_modulus")),
  _damaged_modulus(declareProperty<Real>("damaged_modulus")),
  _a0(declareProperty<Real>("a0")),
  _a1(declareProperty<Real>("a1")),
  _a2(declareProperty<Real>("a2")),
  _a3(declareProperty<Real>("a3")),
  _C_g(declareProperty<Real>("C_g")),
  _m1(declareProperty<Real>("m1")),
  _m2(declareProperty<Real>("m2")), 
  _structural_stress_coefficient(declareProperty<Real>("structural_stress_coefficient")),
  _grad_alpha_damagedvar(declareProperty<RealGradient>("gradient_alpha_damagedvar")),
  _grad_alpha_damagedvar_xdir(declareProperty<Real>("gradient_alpha_damagedvar_xdir")),
  _grad_alpha_damagedvar_ydir(declareProperty<Real>("gradient_alpha_damagedvar_ydir")),
  _velgrad_L(declareProperty<RankTwoTensor>("velgrad_L")),
  //--------------------------------------------------------------//
  //input values
  _lambda_o_value(getParam<Real>("lambda_o")),
  _shear_modulus_o_value(getParam<Real>("shear_modulus_o")),
  _xi_0_value(getParam<Real>("xi_0")),
  _xi_d_value(getParam<Real>("xi_d")),
  _chi_value(getParam<Real>("chi")),
  _C_g_value(getParam<Real>("C_g")),
  _m1_value(getParam<Real>("m1")),
  _m2_value(getParam<Real>("m2")),
  _alpha_damagedvar_aux(coupledValue("alpha_damagedvar_aux")),
  _B_damagedvar_aux(coupledValue("B_damagedvar_aux")),
  _structural_stress_coefficient_aux(coupledValue("structural_stress_coefficient")),
  _grad_alpha_damagedvar_value(coupledGradient("alpha_damagedvar_aux")),
  //--------------------------------------------------------------//
  _grad_vel_x(coupledGradient("vel_x")),
  _grad_vel_y(coupledGradient("vel_y")),
  _grad_vel_z(coupledGradient("vel_z"))
{
}

//Rules:See https://github.com/idaholab/moose/discussions/19450
//Only the object that declares the material property can assign values to it.
//Objects can request material properties, gaining read-only access to their values.
//When any object (including the object that declares it) requests the old value of a material property, that property becomes "stateful".
//All stateful material properties must be initialized within the initQpStatefulProperties call. 
//
void 
DiffusedDamageBreakageMaterialMainApp::initQpStatefulProperties()
{
  /* compute _gamma_damaged_r_mat */
  computegammar();

  /* update damage variable and breakage variable */
  updatedamagebreakage();

  /* compute modulus: _lambda, _shear_modulus, _damaged_modulus */
  updatemodulus();

  /* compute coefficients: a0 a1 a2 a3 */
  computecoefficients();

  /* compute L matrix */
  buildLmatrix();

  /* compute spatial gradient of damage variable */
  _structural_stress_coefficient[_qp] = _structural_stress_coefficient_aux[_qp];
  _grad_alpha_damagedvar[_qp] = _grad_alpha_damagedvar_value[_qp];

  /* compute constant material properties: Cg, m1, m2 */
  _C_g[_qp] = _C_g_value;
  _m1[_qp] = _m1_value;
  _m2[_qp] = _m2_value;
}

void
DiffusedDamageBreakageMaterialMainApp::computeQpProperties()
{
  /* compute _gamma_damaged_r_mat */
  computegammar();

  /* update damage variable and breakage variable */
  updatedamagebreakage();

  /* compute modulus: _lambda, _shear_modulus, _damaged_modulus */
  updatemodulus();

  /* compute coefficients: a0 a1 a2 a3 */
  computecoefficients();

  /* compute L matrix */
  buildLmatrix();

  /* compute spatial gradient of damage variable */
  _structural_stress_coefficient[_qp] = _structural_stress_coefficient_aux[_qp];
  _grad_alpha_damagedvar[_qp] = _grad_alpha_damagedvar_value[_qp];

  //for debugging
  _grad_alpha_damagedvar_xdir[_qp] = _grad_alpha_damagedvar_value[_qp](0);
  _grad_alpha_damagedvar_ydir[_qp] = _grad_alpha_damagedvar_value[_qp](1);

  /* compute constant material properties: Cg, m1, m2 */
  _C_g[_qp] = _C_g_value;
  _m1[_qp] = _m1_value;
  _m2[_qp] = _m2_value;
}

void 
DiffusedDamageBreakageMaterialMainApp::computegammar()
{

  // Calculate each part of the expression
  Real term1 = -_xi_0_value * (-_lambda_o_value * pow(_xi_0_value, 2) + 6 * _lambda_o_value + 2 * _shear_modulus_o_value);
  Real term2_sqrt = sqrt((_lambda_o_value * pow(_xi_0_value, 2) + 2 * _shear_modulus_o_value) * 
                            (_lambda_o_value * pow(_xi_0_value, 4) - 12 * _lambda_o_value * pow(_xi_0_value, 2) + 36 * _lambda_o_value
                            - 6 * _shear_modulus_o_value * pow(_xi_0_value, 2) + 24 * _shear_modulus_o_value));
  Real denominator = 2 * (pow(_xi_0_value, 2) - 3);
  
  // Calculate gamma_r
  Real gamma_r = (term1 - term2_sqrt) / denominator;
  
  //save
  _gamma_damaged_r[_qp] = gamma_r;
}

void
DiffusedDamageBreakageMaterialMainApp::updatedamagebreakage()
{
  _alpha_damagedvar[_qp] = _alpha_damagedvar_aux[_qp];
  _B_damagedvar[_qp] = _B_damagedvar_aux[_qp];
}

void 
DiffusedDamageBreakageMaterialMainApp::updatemodulus()
{
  Real shear_modulus = _shear_modulus_o_value +  _alpha_damagedvar[_qp] * _xi_0_value * _gamma_damaged_r[_qp];
  Real gamma_damaged =  _alpha_damagedvar[_qp] * _gamma_damaged_r[_qp];

  _lambda[_qp] = _lambda_o_value;
  _shear_modulus[_qp] = shear_modulus;
  _damaged_modulus[_qp] = gamma_damaged;
}

void
DiffusedDamageBreakageMaterialMainApp::computecoefficients()
{
  //compute xi_1
  Real _xi_1 = _xi_0_value + sqrt( pow(_xi_0_value , 2) + 2 * _shear_modulus_o_value / _lambda_o_value );

  //compute alpha_cr | xi = 0
  Real alpha_cr_xi0 = alphacr_root1(0);

  //compute mu_cr
  Real mu_cr = _shear_modulus_o_value + alpha_cr_xi0 * _xi_0_value * _gamma_damaged_r[_qp];

  //a0
  Real a0 = _chi_value * mu_cr;

  //a1
  Real numerator_a1 = -2 * _chi_value * mu_cr * pow(_xi_1, 3) + 6 * _chi_value * mu_cr * _xi_1 * pow(_xi_d_value, 2) - 4 * _chi_value * mu_cr * pow(_xi_d_value, 3)
                      - 2 * _gamma_damaged_r[_qp] * pow(_xi_1, 3) * _xi_d_value + 2 * _gamma_damaged_r[_qp] * pow(_xi_1, 3) * _xi_0_value
                      + _lambda_o_value * pow(_xi_1, 3) * pow(_xi_d_value, 2) + 2 * _shear_modulus_o_value * pow(_xi_1, 3);
  Real denominator_a1 = 2 * pow(_xi_1, 3) * _xi_d_value - 4 * pow(_xi_1, 2) * pow(_xi_d_value, 2) + 2 * _xi_1 * pow(_xi_d_value, 3);
  Real a1 = numerator_a1 / denominator_a1;

  //a2
  Real numerator_a2 = 2 * _chi_value * mu_cr * pow(_xi_1, 3) - 3 * _chi_value * mu_cr * pow(_xi_1, 2) * _xi_d_value + _chi_value * mu_cr * pow(_xi_d_value, 3)
                       + 2 * _gamma_damaged_r[_qp] * pow(_xi_1, 3) * _xi_d_value - 2 * _gamma_damaged_r[_qp] * pow(_xi_1, 3) * _xi_0_value
                       - _lambda_o_value * pow(_xi_1, 3) * pow(_xi_d_value, 2) - 2 * _shear_modulus_o_value * pow(_xi_1, 3);
  Real denominator_a2 = pow(_xi_1, 4) * _xi_d_value - 2 * pow(_xi_1, 3) * pow(_xi_d_value, 2) + pow(_xi_1, 2) * pow(_xi_d_value, 3); 
  Real a2 = numerator_a2 / denominator_a2; 

  //a3
  Real numerator_a3 = -2 * _chi_value * mu_cr * pow(_xi_1, 2) + 4 * _chi_value * mu_cr * _xi_1 * _xi_d_value - 2 * _chi_value * mu_cr * pow(_xi_d_value, 2)
                       - 2 * _gamma_damaged_r[_qp] * pow(_xi_1, 2) * _xi_d_value + 2 * _gamma_damaged_r[_qp] * pow(_xi_1, 2) * _xi_0_value
                       + _lambda_o_value * pow(_xi_1, 2) * pow(_xi_d_value, 2) + 2 * _shear_modulus_o_value * pow(_xi_1, 2);
  Real denominator_a3 = 2 * pow(_xi_1, 4) * _xi_d_value - 4 * pow(_xi_1, 3) * pow(_xi_d_value, 2) + 2 * pow(_xi_1, 2) * pow(_xi_d_value, 3);
  Real a3 = numerator_a3 / denominator_a3; 

  //save
  _a0[_qp] = a0;
  _a1[_qp] = a1;
  _a2[_qp] = a2;
  _a3[_qp] = a3;

}

// Function for alpha_func_root1
Real 
DiffusedDamageBreakageMaterialMainApp::alphacr_root1(Real xi) {

  Real term1 = _lambda_o_value * pow(xi, 3) - 6 * _lambda_o_value * _xi_0_value + 6 * _shear_modulus_o_value * xi - 8 * _shear_modulus_o_value * _xi_0_value;
  Real term2 = std::sqrt(_lambda_o_value * _lambda_o_value * pow(xi, 6) 
                            - 12 * _lambda_o_value * _lambda_o_value * pow(xi, 3) * _xi_0_value 
                            + 36 * _lambda_o_value * _lambda_o_value * _xi_0_value * _xi_0_value 
                            + 12 * _lambda_o_value * _shear_modulus_o_value * pow(xi, 4) 
                            - 16 * _lambda_o_value * _shear_modulus_o_value * pow(xi, 3) * _xi_0_value 
                            - 72 * _lambda_o_value * _shear_modulus_o_value * pow(xi, 2) 
                            + 72 * _lambda_o_value * _shear_modulus_o_value * xi * _xi_0_value 
                            + 72 * _lambda_o_value * _shear_modulus_o_value 
                            - 12 * _shear_modulus_o_value * _shear_modulus_o_value * pow(xi, 2) 
                            + 48 * _shear_modulus_o_value * _shear_modulus_o_value);
  Real denominator = 2 * _gamma_damaged_r[_qp] * (3 * pow(xi, 2) - 6 * xi * _xi_0_value + 4 * _xi_0_value * _xi_0_value - 3);
  return (term1 - term2) / denominator;
}

//build the L matrix
//L matrix is the gradient of velocity
void
DiffusedDamageBreakageMaterialMainApp::buildLmatrix()
{
  _velgrad_L[_qp](0,0) = (_grad_vel_x)[_qp](0); _velgrad_L[_qp](0,1) = (_grad_vel_x)[_qp](1); _velgrad_L[_qp](0,2) = (_grad_vel_x)[_qp](2);
  _velgrad_L[_qp](1,0) = (_grad_vel_y)[_qp](0); _velgrad_L[_qp](1,1) = (_grad_vel_y)[_qp](1); _velgrad_L[_qp](1,2) = (_grad_vel_y)[_qp](2);
  _velgrad_L[_qp](2,0) = (_grad_vel_z)[_qp](0); _velgrad_L[_qp](2,1) = (_grad_vel_z)[_qp](1); _velgrad_L[_qp](2,2) = (_grad_vel_z)[_qp](2);
}