//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DamageBreakageMaterial.h"

/**
 *  Material used in damage-breakage large deformation formulation
 *  Created by Chunhui Zhao, Aug 5th, 2024
 */
registerMooseObject("farmsApp", DamageBreakageMaterial);

InputParameters
DamageBreakageMaterial::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Material used in three field poro dynamics simulations");
  params.addRequiredParam<Real>(        "lambda_o", "initial lambda constant value");
  params.addRequiredParam<Real>( "shear_modulus_o", "initial shear modulus value");
  params.addRequiredParam<Real>(            "xi_0", "strain invariants ratio: onset of damage evolution");
  params.addRequiredParam<Real>(            "xi_d", "strain invariants ratio: onset of breakage healing");
  params.addRequiredParam<Real>(          "xi_min", "strain invariants ratio: minimum allowable value");
  params.addRequiredParam<Real>(          "xi_max", "strain invariants ratio: maximum allowable value");
  params.addRequiredParam<Real>(             "chi", "coefficient of energy ratio Fb/Fs = chi < 1");
  params.addRequiredParam<Real>(     "Cd_constant", "coefficient gives positive damage evolution");
  params.addRequiredParam<Real>(             "C_1", "coefficient of healing for damage evolution");
  params.addRequiredParam<Real>(             "C_2", "coefficient of healing for damage evolution");
  params.addRequiredParam<Real>(      "beta_width", "coefficient gives width of transitional region");
  params.addRequiredParam<Real>( "CdCb_multiplier", "multiplier between Cd and Cb");
  params.addRequiredParam<Real>(    "CBH_constant", "constant CBH value");
  params.addRequiredParam<Real>(             "C_g", "compliance or fluidity of the fine grain granular material");
  params.addRequiredParam<Real>(              "m1", "coefficient of power law indexes");
  params.addRequiredParam<Real>(              "m2", "coefficient of power law indexes");
  params.addParam<bool>("use_xi0_aux", false, "Whether to use xi_0 from aux variable");
  params.addCoupledVar("xi0_aux", "AuxVariable for xi_0");
  return params;
}

DamageBreakageMaterial::DamageBreakageMaterial(const InputParameters & parameters)
  : Material(parameters),
  _alpha_damagedvar(declareProperty<Real>("alpha_damagedvar")),
  _B_breakagevar(declareProperty<Real>("B_damagedvar")),
  _lambda(declareProperty<Real>("lambda_const")),
  _shear_modulus(declareProperty<Real>("shear_modulus")),
  _damaged_modulus(declareProperty<Real>("damaged_modulus")),
  _gamma_damaged_r(declareProperty<Real>("gamma_damaged_r")),
  _a0(declareProperty<Real>("a0")),
  _a1(declareProperty<Real>("a1")),
  _a2(declareProperty<Real>("a2")),
  _a3(declareProperty<Real>("a3")),
  _C_g(declareProperty<Real>("C_g")),
  _m1(declareProperty<Real>("m1")),
  _m2(declareProperty<Real>("m2")),  
  _xi_1_mat(declareProperty<Real>("xi_1")),  
  _xi_0_mat(declareProperty<Real>("xi_0")),
  _alpha_damagedvar_old(getMaterialPropertyOldByName<Real>("alpha_damagedvar")),
  _B_breakagevar_old(getMaterialPropertyOldByName<Real>("B_damagedvar")),
  _I2_old(getMaterialPropertyOldByName<Real>("second_elastic_strain_invariant")),
  _xi_old(getMaterialPropertyOldByName<Real>("strain_invariant_ratio")),
  _initial_damage(getMaterialProperty<Real>("initial_damage")),
  _gamma_damaged_r_old(getMaterialPropertyOldByName<Real>("gamma_damaged_r")),
  _a0_old(getMaterialPropertyOldByName<Real>("a0")),
  _a1_old(getMaterialPropertyOldByName<Real>("a1")),
  _a2_old(getMaterialPropertyOldByName<Real>("a2")),
  _a3_old(getMaterialPropertyOldByName<Real>("a3")),
  _lambda_o(getParam<Real>("lambda_o")),
  _shear_modulus_o(getParam<Real>("shear_modulus_o")),
  _xi_d(getParam<Real>("xi_d")),
  _xi_min(getParam<Real>("xi_min")),
  _xi_max(getParam<Real>("xi_max")),
  _chi(getParam<Real>("chi")),
  _Cd_constant(getParam<Real>("Cd_constant")), 
  _C1(getParam<Real>("C_1")),
  _C2(getParam<Real>("C_2")),
  _beta_width(getParam<Real>("beta_width")),
  _CdCb_multiplier(getParam<Real>("CdCb_multiplier")),
  _CBH_constant(getParam<Real>("CBH_constant")),
  _C_g_value(getParam<Real>("C_g")),
  _m1_value(getParam<Real>("m1")),
  _m2_value(getParam<Real>("m2")),
  _use_xi0_aux(getParam<bool>("use_xi0_aux")),
  _xi0_value(getParam<Real>("xi_0")),
  _xi0_aux(_use_xi0_aux ? &coupledValue("xi0_aux") : nullptr)
{
  if (_use_xi0_aux && !parameters.isParamSetByUser("xi0_aux"))
    mooseError("Must specify xi0_aux when use_xi0_aux = true");
}

//Rules:See https://github.com/idaholab/moose/discussions/19450
//Only the object that declares the material property can assign values to it.
//Objects can request material properties, gaining read-only access to their values.
//When any object (including the object that declares it) requests the old value of a material property, that property becomes "stateful".
//All stateful material properties must be initialized within the initQpStatefulProperties call. 
//
void 
DamageBreakageMaterial::initQpStatefulProperties()
{

  // Get xi_0 value
  const Real _xi_0 = _use_xi0_aux ? (*_xi0_aux)[_qp] : _xi0_value;

  /* compute gamma_r */
  computegammar();

  /* compute a0 a1 a2 a3 coefficients */
  computecoefficients();

  _alpha_damagedvar[_qp] = _initial_damage[_qp];
  _B_breakagevar[_qp] = 0.0;

  _lambda[_qp] = _lambda_o;
  _shear_modulus[_qp] = _shear_modulus_o + _initial_damage[_qp] * _xi_0 * _gamma_damaged_r[_qp];
  _damaged_modulus[_qp] = _initial_damage[_qp] * _gamma_damaged_r[_qp];
  
  /* define Cg m1 m2 */
  _m1[_qp] = _m1_value;
  _m2[_qp] = _m2_value;
  _C_g[_qp] = _C_g_value;

}

void
DamageBreakageMaterial::computeQpProperties()
{

  /* compute gamma_r */
  computegammar();

  /* compute a0 a1 a2 a3 coefficients */
  computecoefficients();

  /* compute alpha at t_{n+1} using quantities from t_{n} */
  Real alpha_updated = updatedamage();

  /* compute B at t_{n+1} using quantities from t_{n} */
  updatebreakage();

  /* compute modulus at t_{n+1} using alpha at t_{n} */
  updatemodulus(alpha_updated);

  /* define Cg m1 m2 */
  _m1[_qp] = _m1_value;
  _m2[_qp] = _m2_value;
  _C_g[_qp] = _C_g_value;

  // Get xi_0 value
  const Real _xi_0 = _use_xi0_aux ? (*_xi0_aux)[_qp] : _xi0_value;
  _xi_0_mat[_qp] = _xi_0;

}

Real 
DamageBreakageMaterial::updatedamage()
{

  // Get xi_0 value
  const Real _xi_0 = _use_xi0_aux ? (*_xi0_aux)[_qp] : _xi0_value;

  //compute forcing term
  Real alpha_forcingterm;
  if ( _xi_old[_qp] >= _xi_0 && _xi_old[_qp] <= _xi_max ){
    alpha_forcingterm = (1 - _B_breakagevar_old[_qp]) * ( _Cd_constant * _I2_old[_qp] * ( _xi_old[_qp] - _xi_0 ) );
  }
  else if ( _xi_old[_qp] < _xi_0 && _xi_old[_qp] >= _xi_min ){
    alpha_forcingterm = (1 - _B_breakagevar_old[_qp]) * ( _C1 * std::exp(_alpha_damagedvar_old[_qp]/_C2) * _I2_old[_qp] * ( _xi_old[_qp] - _xi_0 ) );
  }
  else{
    //mooseError("xi_old is OUT-OF-RANGE!.");   
  }

  //update alpha at current time
  Real alpha_damagedvar = _alpha_damagedvar_old[_qp] + _dt * alpha_forcingterm;

  //check alpha within range
  if ( alpha_damagedvar < 0 ){ alpha_damagedvar = 0.0; }
  else if ( alpha_damagedvar > 1 ){ alpha_damagedvar = 1.0; }
  else{} 

  //check below initial damage (fix initial damage)
  if ( alpha_damagedvar < _initial_damage[_qp] ){ alpha_damagedvar = _initial_damage[_qp]; }
  else{}

  _alpha_damagedvar[_qp] = alpha_damagedvar;

  return alpha_damagedvar;
}

void
DamageBreakageMaterial::updatebreakage()
{

  // Get xi_0 value
  const Real _xi_0 = _use_xi0_aux ? (*_xi0_aux)[_qp] : _xi0_value;

  /* compute C_B based on C_d */
  Real C_B = _CdCb_multiplier * _Cd_constant;

  //compute xi_1
  Real _xi_1 = _xi_0 + sqrt( pow(_xi_0 , 2) + 2 * _shear_modulus_o / _lambda_o );

  //save
  _xi_1_mat[_qp] = _xi_1;

  //alphacr function
  Real alphacr;
  if ( _xi_old[_qp] < _xi_0 ){ alphacr = 1.0;} 
  else if ( _xi_old[_qp] > _xi_0 && _xi_old[_qp] <= _xi_1 ){ alphacr = alphacr_root1(_xi_old[_qp]);}
  else if ( _xi_old[_qp] > _xi_1 && _xi_old[_qp] <= _xi_max ){ alphacr = alphacr_root2(_xi_old[_qp]);}
  else{std::cout<<"xi: "<<_xi_old[_qp]<<std::endl;}//mooseError("xi exceeds the maximum allowable range!")}

  //compute forcing func
  Real Prob = 1.0 / ( std::exp( (alphacr - _alpha_damagedvar_old[_qp]) / _beta_width ) + 1.0 );
  Real B_forcingterm;
  if ( _xi_old[_qp] >= _xi_d && _xi_old[_qp] <= _xi_max ){
    B_forcingterm = 1.0 * C_B * Prob * (1-_B_breakagevar_old[_qp]) * _I2_old[_qp] * (_xi_old[_qp] - _xi_d); //could heal if xi < xi_0
  }
  else if ( _xi_old[_qp] < _xi_d && _xi_old[_qp] >= _xi_min ){
    B_forcingterm = 1.0 * _CBH_constant * _I2_old[_qp] * ( _xi_old[_qp] - _xi_d );
  }
  else{
    //mooseError("xi_old is OUT-OF-RANGE!.");
  }

  Real B_damagedvar = _B_breakagevar_old[_qp] + _dt * B_forcingterm;

  //check breakage within range
  if ( B_damagedvar < 0 ){ B_damagedvar = 0.0; }
  else if ( B_damagedvar > 1 ){ B_damagedvar = 1.0; }
  else{}   

  _B_breakagevar[_qp] = B_damagedvar;

}

void 
DamageBreakageMaterial::updatemodulus(Real alpha_updated)
{

  // Get xi_0 value
  const Real _xi_0 = _use_xi0_aux ? (*_xi0_aux)[_qp] : _xi0_value;

  Real shear_modulus = _shear_modulus_o + alpha_updated * _xi_0 * _gamma_damaged_r[_qp];
  Real gamma_damaged = alpha_updated * _gamma_damaged_r[_qp];

  _lambda[_qp] = _lambda_o;
  _shear_modulus[_qp] = shear_modulus;
  _damaged_modulus[_qp] = gamma_damaged;
}

void 
DamageBreakageMaterial::computegammar()
{

  // Get xi_0 value
  const Real _xi_0 = _use_xi0_aux ? (*_xi0_aux)[_qp] : _xi0_value;

  // Calculate each part of the expression
  Real term1 = -_xi_0 * (-_lambda_o * pow(_xi_0, 2) + 6 * _lambda_o + 2 * _shear_modulus_o);
  Real term2_sqrt = sqrt((_lambda_o * pow(_xi_0, 2) + 2 * _shear_modulus_o) * 
                            (_lambda_o * pow(_xi_0, 4) - 12 * _lambda_o * pow(_xi_0, 2) + 36 * _lambda_o 
                            - 6 * _shear_modulus_o * pow(_xi_0, 2) + 24 * _shear_modulus_o));
  Real denominator = 2 * (pow(_xi_0, 2) - 3);
  
  // Calculate gamma_r
  Real gamma_r = (term1 - term2_sqrt) / denominator;
  
  //save
  _gamma_damaged_r[_qp] = gamma_r;
}

void
DamageBreakageMaterial::computecoefficients()
{

  // Get xi_0 value
  const Real _xi_0 = _use_xi0_aux ? (*_xi0_aux)[_qp] : _xi0_value;

  //compute xi_1
  Real _xi_1 = _xi_0 + sqrt( pow(_xi_0 , 2) + 2 * _shear_modulus_o / _lambda_o );

  //compute alpha_cr | xi = 0
  Real alpha_cr_xi0 = alphacr_root1(0);

  //compute mu_cr
  Real mu_cr = _shear_modulus_o + alpha_cr_xi0 * _xi_0 * _gamma_damaged_r[_qp];

  //a0
  Real a0 = _chi * mu_cr;

  //a1
  Real numerator_a1 = -2 * _chi * mu_cr * pow(_xi_1, 3) + 6 * _chi * mu_cr * _xi_1 * pow(_xi_d, 2) - 4 * _chi * mu_cr * pow(_xi_d, 3)
                      - 2 * _gamma_damaged_r[_qp] * pow(_xi_1, 3) * _xi_d + 2 * _gamma_damaged_r[_qp] * pow(_xi_1, 3) * _xi_0
                      + _lambda_o * pow(_xi_1, 3) * pow(_xi_d, 2) + 2 * _shear_modulus_o * pow(_xi_1, 3);
  Real denominator_a1 = 2 * pow(_xi_1, 3) * _xi_d - 4 * pow(_xi_1, 2) * pow(_xi_d, 2) + 2 * _xi_1 * pow(_xi_d, 3);
  Real a1 = numerator_a1 / denominator_a1;

  //a2
  Real numerator_a2 = 2 * _chi * mu_cr * pow(_xi_1, 3) - 3 * _chi * mu_cr * pow(_xi_1, 2) * _xi_d + _chi * mu_cr * pow(_xi_d, 3)
                       + 2 * _gamma_damaged_r[_qp] * pow(_xi_1, 3) * _xi_d - 2 * _gamma_damaged_r[_qp] * pow(_xi_1, 3) * _xi_0
                       - _lambda_o * pow(_xi_1, 3) * pow(_xi_d, 2) - 2 * _shear_modulus_o * pow(_xi_1, 3);
  Real denominator_a2 = pow(_xi_1, 4) * _xi_d - 2 * pow(_xi_1, 3) * pow(_xi_d, 2) + pow(_xi_1, 2) * pow(_xi_d, 3); 
  Real a2 = numerator_a2 / denominator_a2; 

  //a3
  Real numerator_a3 = -2 * _chi * mu_cr * pow(_xi_1, 2) + 4 * _chi * mu_cr * _xi_1 * _xi_d - 2 * _chi * mu_cr * pow(_xi_d, 2)
                       - 2 * _gamma_damaged_r[_qp] * pow(_xi_1, 2) * _xi_d + 2 * _gamma_damaged_r[_qp] * pow(_xi_1, 2) * _xi_0
                       + _lambda_o * pow(_xi_1, 2) * pow(_xi_d, 2) + 2 * _shear_modulus_o * pow(_xi_1, 2);
  Real denominator_a3 = 2 * pow(_xi_1, 4) * _xi_d - 4 * pow(_xi_1, 3) * pow(_xi_d, 2) + 2 * pow(_xi_1, 2) * pow(_xi_d, 3);
  Real a3 = numerator_a3 / denominator_a3; 

  //save
  _a0[_qp] = a0;
  _a1[_qp] = a1;
  _a2[_qp] = a2;
  _a3[_qp] = a3;

}

// Function for alpha_func_root1
Real 
DamageBreakageMaterial::alphacr_root1(Real xi) {
    
    // Get xi_0 value
    const Real _xi_0 = _use_xi0_aux ? (*_xi0_aux)[_qp] : _xi0_value;

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
    Real denominator = 2 * _gamma_damaged_r[_qp] * (3 * pow(xi, 2) - 6 * xi * _xi_0 + 4 * _xi_0 * _xi_0 - 3);
    return (term1 - term2) / denominator;
}

// Function for alpha_func_root2
Real 
DamageBreakageMaterial::alphacr_root2(Real xi) {
    
    // Get xi_0 value
    const Real _xi_0 = _use_xi0_aux ? (*_xi0_aux)[_qp] : _xi0_value;

    return 2 * _shear_modulus_o / (_gamma_damaged_r[_qp] * (xi - 2 * _xi_0));
}