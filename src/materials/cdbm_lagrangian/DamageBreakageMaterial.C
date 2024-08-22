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
  params.addRequiredParam<Real>(            "xi_1", "critical point of three phases");
  params.addRequiredParam<Real>(          "xi_min", "strain invariants ratio: minimum allowable value");
  params.addRequiredParam<Real>(          "xi_max", "strain invariants ratio: maximum allowable value");
  params.addRequiredParam<Real>( "gamma_damaged_r", "coefficient of damage solid modulus");
  params.addRequiredParam<Real>(              "m1", "coefficient of std::power law indexes");
  params.addRequiredParam<Real>(              "m2", "coefficient of std::power law indexes");
  params.addRequiredParam<Real>(     "Cd_constant", "coefficient gives positive damage evolution");
  params.addRequiredParam<Real>(             "C_1", "coefficient of healing for damage evolution");
  params.addRequiredParam<Real>(             "C_2", "coefficient of healing for damage evolution");
  params.addRequiredParam<Real>(      "beta_width", "coefficient gives width of transitional region");
  params.addRequiredParam<Real>( "CdCb_multiplier", "multiplier between Cd and Cb");
  params.addRequiredParam<Real>(    "CBH_constant", "constant CBH value");
  return params;
}

DamageBreakageMaterial::DamageBreakageMaterial(const InputParameters & parameters)
  : Material(parameters),
  _alpha_damagedvar(declareProperty<Real>("alpha_damagedvar")),
  _B_breakagevar(declareProperty<Real>("B_damagedvar")),
  _shear_modulus(declareProperty<Real>("shear_modulus")),
  _damaged_modulus(declareProperty<Real>("damaged_modulus")),
  _alpha_damagedvar_old(getMaterialPropertyOldByName<Real>("alpha_damagedvar")),
  _B_breakagevar_old(getMaterialPropertyOldByName<Real>("B_damagedvar")),
  _I2_old(getMaterialPropertyOldByName<Real>("I2")),
  _xi_old(getMaterialPropertyOldByName<Real>("xi")),
  _initial_damage(getMaterialProperty<Real>("initial_damage")),
  _shear_modulus_o(getParam<Real>("shear_modulus_o")),
  _xi_0(getParam<Real>("xi_0")),
  _xi_d(getParam<Real>("xi_d")),
  _xi_1(getParam<Real>("xi_1")),
  _xi_min(getParam<Real>("xi_min")),
  _xi_max(getParam<Real>("xi_max")),
  _gamma_damaged_r(getParam<Real>("gamma_damaged_r")),
  _Cd_constant(getParam<Real>("Cd_constant")), 
  _C1(getParam<Real>("C_1")),
  _C2(getParam<Real>("C_2")),
  _beta_width(getParam<Real>("beta_width")),
  _CdCb_multiplier(getParam<Real>("CdCb_multiplier")),
  _CBH_constant(getParam<Real>("CBH_constant"))
{
}

void
DamageBreakageMaterial::computeQpProperties()
{
  /* compute alpha at t_{n+1} using quantities from t_{n} */
  Real alpha_updated = updatedamage();

  /* compute B at t_{n+1} using quantities from t_{n} */
  updatebreakage();

  /* compute modulus at t_{n+1} using alpha at t_{n} */
  updatemodulus(alpha_updated);

}

Real 
DamageBreakageMaterial::updatedamage()
{
  //compute forcing term
  Real alpha_forcingterm;
  if ( _xi_old[_qp] >= _xi_0 && _xi_old[_qp] <= _xi_max ){
    alpha_forcingterm = (1 - _B_breakagevar_old[_qp]) * ( _Cd_constant * _I2_old[_qp] * ( _xi_old[_qp] - _xi_0 ) );
  }
  else if ( _xi_old[_qp] < _xi_0 && _xi_old[_qp] >= _xi_min ){
    alpha_forcingterm = (1 - _B_breakagevar_old[_qp]) * ( _C1 * std::exp(_alpha_damagedvar_old[_qp]/_C2) * _I2_old[_qp] * ( _xi_old[_qp] - _xi_0 ) );
  }
  else{
    mooseError("xi_old is OUT-OF-RANGE!.");   
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
  /* compute C_B based on C_d */
  Real C_B = _CdCb_multiplier * _Cd_constant;

  //alphacr function
  Real alphacr;
  if ( _xi_old[_qp] < _xi_0 ){ alphacr = 1.0;} 
  else if ( _xi_old[_qp] > _xi_0 && _xi_old[_qp] <= _xi_1 ){ alphacr = ((_xi_old[_qp]*2.76e5-7.100521107637101e2*_xi_old[_qp]*7.5e2-7.100521107637101e2*1.4e3-7.100521107637101e2*std::pow(_xi_old[_qp],3)*1.25e2+std::pow(_xi_old[_qp],3)*4.6e4+std::sqrt((7.100521107637101e2*3.68e2-3.19799e5)*(_xi_old[_qp]*-1.44e3-std::pow(_xi_old[_qp],2)*2.1e3+std::pow(_xi_old[_qp],3)*5.6e2+std::pow(_xi_old[_qp],4)*3.0e2+std::pow(_xi_old[_qp],6)*2.5e1+3.576e3)*(-3.590922148807814e-1))*5.9e1+5.152e5)*(-5.9e1/4.0))/(_xi_old[_qp]*3.837588e7-7.100521107637101e2*_xi_old[_qp]*4.416e4+7.100521107637101e2*4.048e3-7.100521107637101e2*std::pow(_xi_old[_qp],2)*2.76e4+std::pow(_xi_old[_qp],2)*2.3984925e7-3.517789e6);}
  else if ( _xi_old[_qp] > _xi_1 && _xi_old[_qp] <= _xi_max ){ alphacr = 6.0e10/(7.100521107637101e2*1.627118644067797e8+_xi_old[_qp]*(7.100521107637101e2*1.016949152542373e8-3.742372881355932e10)-5.987796610169492e10);}
  else{std::cout<<"xi: "<<_xi_old[_qp]<<std::endl;mooseError("xi exceeds the maximum allowable range!");}

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
    mooseError("xi_old is OUT-OF-RANGE!.");
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
  Real shear_modulus = _shear_modulus_o + alpha_updated * _xi_0 * _gamma_damaged_r;
  Real gamma_damaged = alpha_updated * _gamma_damaged_r;
  _shear_modulus[_qp] = shear_modulus;
  _damaged_modulus[_qp] = gamma_damaged;
}