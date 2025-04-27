//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DiffusedDamageBreakageMaterialSubApp.h"

/**
 *  Material used in damage-breakage large deformation formulation, consider full damage evolution equation with diffusion
 *  Created by Chunhui Zhao, Dec 24th, 2024
 */
registerMooseObject("farmsApp", DiffusedDamageBreakageMaterialSubApp);

InputParameters
DiffusedDamageBreakageMaterialSubApp::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Material used in three field poro dynamics simulations");
  //input parameters
  params.addRequiredParam<Real>(        "lambda_o", "initial lambda constant value");
  params.addRequiredParam<Real>( "shear_modulus_o", "initial shear modulus value");
  params.addRequiredParam<Real>(            "xi_0", "strain invariants ratio: onset of damage evolution");
  params.addRequiredParam<Real>(            "xi_d", "strain invariants ratio: onset of breakage healing");
  params.addRequiredParam<Real>(          "xi_min", "strain invariants ratio: minimum allowable value");
  params.addRequiredParam<Real>(          "xi_max", "strain invariants ratio: maximum allowable value");
  params.addRequiredParam<Real>(     "Cd_constant", "coefficient gives positive damage evolution");
  params.addRequiredParam<Real>(             "C_1", "coefficient of healing for damage evolution");
  params.addRequiredParam<Real>(             "C_2", "coefficient of healing for damage evolution");
  params.addRequiredParam<Real>(      "beta_width", "coefficient gives width of transitional region");
  params.addRequiredParam<Real>( "CdCb_multiplier", "multiplier between Cd and Cb");
  params.addRequiredParam<Real>(    "CBH_constant", "constant CBH value");
  params.addRequiredParam<Real>(     "D_diffusion", "material parameter: compliance or fluidity of the fine grain granular material");
  //input coupled variables from main app
  params.addRequiredCoupledVar("I2_aux", "second_elastic_strain_invariant");
  params.addRequiredCoupledVar("xi_aux", "strain_invariant_ratio");
  params.addRequiredCoupledVar("initial_damage_aux", "initial_damage");
  //get strain rate dependent Cd options
  params.addParam<bool>("use_cd_strain_dependent", false, "Use strain rate dependent Cd");
  params.addParam<Real>("m_exponent", -1.0, "Exponent for strain rate dependent Cd");
  params.addParam<Real>("strain_rate_hat", -1.0, "Strain rate hat for strain rate dependent Cd");
  params.addParam<Real>("cd_hat", -1.0, "Cd hat for strain rate dependent Cd");
  params.addCoupledVar("strain_rate", "strain_rate");
  return params;
}

DiffusedDamageBreakageMaterialSubApp::DiffusedDamageBreakageMaterialSubApp(const InputParameters & parameters)
  : Material(parameters),
  //declare properties
  //--------------------------------------------------------------//
  _lambda_o_mat(declareProperty<Real>("lambda_o")),
  _shear_modulus_o_mat(declareProperty<Real>("shear_modulus_o")),
  _gamma_damaged_r_mat(declareProperty<Real>("gamma_damaged_r")), 
  _xi_1_mat(declareProperty<Real>("xi_1")),  
  _xi_0_mat(declareProperty<Real>("xi_0")),
  _xi_d_mat(declareProperty<Real>("xi_d")),
  _xi_min_mat(declareProperty<Real>("xi_min")),
  _xi_max_mat(declareProperty<Real>("xi_max")),
  _Cd_mat(declareProperty<Real>("Cd")),
  _CdCb_multiplier_mat(declareProperty<Real>("CdCb_multiplier")),
  _CBH_constant_mat(declareProperty<Real>("CBH_constant")),
  _beta_width_mat(declareProperty<Real>("beta_width")),
  _C1_mat(declareProperty<Real>("C_1")),
  _C2_mat(declareProperty<Real>("C_2")),
  _D_diffusion_mat(declareProperty<Real>("D_diffusion")),
  _initial_damage_mat(declareProperty<Real>("initial_damage")),
  _I2_mat(declareProperty<Real>("I2")),
  _xi_mat(declareProperty<Real>("xi")),
  _structural_stress_coefficient_mat(declareProperty<Real>("structural_stress_coefficient")),
  //--------------------------------------------------------------//
  //input values
  _lambda_o_value(getParam<Real>("lambda_o")),
  _shear_modulus_o_value(getParam<Real>("shear_modulus_o")),
  _xi_0_value(getParam<Real>("xi_0")),
  _xi_d_value(getParam<Real>("xi_d")),
  _xi_min_value(getParam<Real>("xi_min")),
  _xi_max_value(getParam<Real>("xi_max")),
  _Cd_constant_value(getParam<Real>("Cd_constant")),
  _CdCb_multiplier_value(getParam<Real>("CdCb_multiplier")),
  _CBH_constant_value(getParam<Real>("CBH_constant")), 
  _beta_width_value(getParam<Real>("beta_width")),
  _C1_value(getParam<Real>("C_1")),
  _C2_value(getParam<Real>("C_2")),
  _D_diffusion_value(getParam<Real>("D_diffusion")),
  _I2_aux(coupledValue("I2_aux")),
  _xi_aux(coupledValue("xi_aux")),
  _initial_damage_aux(coupledValue("initial_damage_aux")),
  //--------------------------------------------------------------//
  //strain rate dependent Cd options
  _use_cd_strain_dependent(getParam<bool>("use_cd_strain_dependent")),
  _strain_rate_hat(getParam<Real>("strain_rate_hat")),
  _cd_hat(getParam<Real>("cd_hat")),
  _m_exponent(getParam<Real>("m_exponent")),
  _strain_rate(coupledValue("strain_rate"))
  //--------------------------------------------------------------//
{
  //check strain rate dependent Cd options
  if (_use_cd_strain_dependent && (_strain_rate_hat < 0 || _cd_hat < 0 || _m_exponent < 0)){
    mooseError("Strain rate dependent Cd options are not set correctly");
  }
  if (!_use_cd_strain_dependent && (_Cd_constant_value < 0)){
    mooseError("Cd_constant is not set correctly, the strain rate dependent is not on");
  }
}

//Rules:See https://github.com/idaholab/moose/discussions/19450
//Only the object that declares the material property can assign values to it.
//Objects can request material properties, gaining read-only access to their values.
//When any object (including the object that declares it) requests the old value of a material property, that property becomes "stateful".
//All stateful material properties must be initialized within the initQpStatefulProperties call. 
//
void 
DiffusedDamageBreakageMaterialSubApp::initQpStatefulProperties()
{
  /* copmute _lambda_o_mat, _shear_modulus_o_mat*/
  _lambda_o_mat[_qp] = _lambda_o_value;
  _shear_modulus_o_mat[_qp] = _shear_modulus_o_value;

  /* compute _gamma_damaged_r_mat */
  computegammar();

  /* compute _xi_1_mat, xi_0_mat, _xi_d_mat, _xi_min_mat, _xi_max_mat */
  _xi_1_mat[_qp] = _xi_0_value + sqrt( pow(_xi_0_value , 2) + 2 * _shear_modulus_o_value / _lambda_o_value );
  _xi_0_mat[_qp] = _xi_0_value;
  _xi_d_mat[_qp] = _xi_d_value;
  _xi_min_mat[_qp] = _xi_min_value;
  _xi_max_mat[_qp] = _xi_max_value;

  /* compute _Cd_mat, _CdCb_multiplier_mat, _CBH_constant_mat */
  _Cd_mat[_qp] = _Cd_constant_value;
  _CdCb_multiplier_mat[_qp] = _CdCb_multiplier_value;
  _CBH_constant_mat[_qp] = _CBH_constant_value;

  /* compute _beta_width_mat, _C1_mat, _C2_mat */
  _beta_width_mat[_qp] = _beta_width_value;
  _C1_mat[_qp] = _C1_value;
  _C2_mat[_qp] = _C2_value;

  /* compute _D_diffusion_mat */
  _D_diffusion_mat[_qp] = _D_diffusion_value;

  /* compute _initial_damage_mat, _I2_mat, _xi_mat */
  _I2_mat[_qp] = _I2_aux[_qp];
  _xi_mat[_qp] = _xi_aux[_qp];
  _initial_damage_mat[_qp] = _initial_damage_aux[_qp];

  /* compute structural stress coefficient */
  _structural_stress_coefficient_mat[_qp] = _D_diffusion_value * _shear_modulus_o_value / _Cd_constant_value;
}

void
DiffusedDamageBreakageMaterialSubApp::computeQpProperties()
{

  /* copmute _lambda_o_mat, _shear_modulus_o_mat*/
  _lambda_o_mat[_qp] = _lambda_o_value;
  _shear_modulus_o_mat[_qp] = _shear_modulus_o_value;

  /* compute _gamma_damaged_r_mat */
  computegammar();

  /* compute _xi_1_mat, xi_0_mat, _xi_d_mat, _xi_min_mat, _xi_max_mat */
  _xi_1_mat[_qp] = _xi_0_value + sqrt( pow(_xi_0_value , 2) + 2 * _shear_modulus_o_value / _lambda_o_value );
  _xi_0_mat[_qp] = _xi_0_value;
  _xi_d_mat[_qp] = _xi_d_value;
  _xi_min_mat[_qp] = _xi_min_value;
  _xi_max_mat[_qp] = _xi_max_value;

  /* compute _Cd_mat, _CdCb_multiplier_mat, _CBH_constant_mat */
  if (_use_cd_strain_dependent){ //option to use strain rate dependent
    computeStrainRateCd();
  }
  else{
    _Cd_mat[_qp] = _Cd_constant_value;
  }
  
  _CdCb_multiplier_mat[_qp] = _CdCb_multiplier_value;
  _CBH_constant_mat[_qp] = _CBH_constant_value;

  /* compute _beta_width_mat, _C1_mat, _C2_mat */
  _beta_width_mat[_qp] = _beta_width_value;
  _C1_mat[_qp] = _C1_value;
  _C2_mat[_qp] = _C2_value;

  /* compute _D_diffusion_mat */
  _D_diffusion_mat[_qp] = _D_diffusion_value;

  /* compute _initial_damage_mat, _I2_mat, _xi_mat */
  _I2_mat[_qp] = _I2_aux[_qp];
  _xi_mat[_qp] = _xi_aux[_qp];
  _initial_damage_mat[_qp] = _initial_damage_aux[_qp];

  /* compute structural stress coefficient */
  // if (_Cd_mat[_qp] == 0){
  //   _structural_stress_coefficient_mat[_qp] = _D_diffusion_value * _shear_modulus_o_value / _cd_hat;
  // }
  // else{
  //   _structural_stress_coefficient_mat[_qp] = _D_diffusion_value * _shear_modulus_o_value / _Cd_mat[_qp];
  // }
  _structural_stress_coefficient_mat[_qp] = 0.0;
}

void 
DiffusedDamageBreakageMaterialSubApp::computegammar()
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
  _gamma_damaged_r_mat[_qp] = gamma_r;
}

void 
DiffusedDamageBreakageMaterialSubApp::computeStrainRateCd()
{
  //_m_exponent: constant value - default value = 0.8
  //_strain_rate_hat: constant value - default value = 1e-4
  //_cd_hat: constant value - default value = 1
  //_strain_rate: deviatoric strain rate, variable value passed from main app
  _Cd_mat[_qp] = pow(10, 1 + _m_exponent * std::log10(_strain_rate[_qp]/_strain_rate_hat)) * _cd_hat;
}