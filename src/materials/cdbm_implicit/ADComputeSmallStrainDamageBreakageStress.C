//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADComputeSmallStrainDamageBreakageStress.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include "SymmetricRankTwoTensor.h"
#include "SymmetricRankFourTensor.h"

#include <random>

registerMooseObject("farmsApp", ADComputeSmallStrainDamageBreakageStress);

InputParameters
ADComputeSmallStrainDamageBreakageStress::validParams()
{
  InputParameters params = ADComputeDamageBreakageStressBase::validParams();
  params.addClassDescription("Compute stress using elasticity for small strains");

  //constant parameters
  params.addRequiredParam<Real>(        "lambda_o", "initial lambda constant value");
  params.addRequiredParam<Real>( "shear_modulus_o", "initial shear modulus value");
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
  params.addRequiredParam<Real>(     "Cd_constant", "coefficient gives positive damage evolution");
  params.addRequiredParam<Real>(             "C_1", "coefficient of healing for damage evolution");
  params.addRequiredParam<Real>(             "C_2", "coefficient of healing for damage evolution");
  params.addRequiredParam<Real>(      "beta_width", "coefficient gives width of transitional region");
  params.addRequiredParam<Real>( "CdCb_multiplier", "multiplier between Cd and Cb");
  params.addRequiredParam<Real>(    "CBH_constant", "constant CBH value");
  return params;
}

ADComputeSmallStrainDamageBreakageStress::ADComputeSmallStrainDamageBreakageStress(
    const InputParameters & parameters)
  : ADComputeDamageBreakageStressBase(parameters),
    _lambda_o(getParam<Real>("lambda_o")),
    _shear_modulus_o(getParam<Real>("shear_modulus_o")),
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
    _Cd_constant(getParam<Real>("Cd_constant")),
    _C1(getParam<Real>("C_1")),
    _C2(getParam<Real>("C_2")),
    _beta_width(getParam<Real>("beta_width")),
    _CdCb_multiplier(getParam<Real>("CdCb_multiplier")),
    _CBH_constant(getParam<Real>("CBH_constant")),
    _initial_damage(getMaterialPropertyOld<Real>("initial_damage")),
    _alpha_damagedvar_old(getMaterialPropertyOld<Real>("alpha_damagedvar")),
    _B_breakagevar_old(getMaterialPropertyOld<Real>("B_breakagevar")),
    _xi_old(getMaterialPropertyOld<Real>("xi")),
    _I2_old(getMaterialPropertyOld<Real>("I2")),
    _shear_modulus_old(getMaterialPropertyOld<Real>("shear_modulus")),
    _gamma_damaged_old(getMaterialPropertyOld<Real>("gamma_damaged")),
    _eps_p_old(getMaterialPropertyOld<RankTwoTensor>("eps_p")),
    _eps_e_old(getMaterialPropertyOld<RankTwoTensor>("eps_e")),
    _sigma_d_old(getMaterialPropertyOld<RankTwoTensor>("sigma_d")),
    _step(_fe_problem.timeStep())
{
}

void
ADComputeSmallStrainDamageBreakageStress::initialSetup()
{
  if ( hasBlockMaterialProperty<RankTwoTensor>(_base_name + "strain_increment"))
    mooseError("This linear elastic stress calculation only works for small strains; use "
               "ADComputeFiniteStrainElasticStress for simulations using incremental and finite "
               "strains.");
}

void
ADComputeSmallStrainDamageBreakageStress::computeQpStress()
{
    if (_step == 1){

    _alpha_damagedvar[_qp] = _initial_damage[_qp];
    _B_breakagevar[_qp] = 0.0;

    /* update modulus */
    ADReal shear_modulus = _shear_modulus_o + _initial_damage[_qp] * _xi_0 * _gamma_damaged_r;
    ADReal gamma_damaged = _initial_damage[_qp] * _gamma_damaged_r;
    _shear_modulus[_qp] = shear_modulus;
    _gamma_damaged[_qp] = gamma_damaged;  

    // stress = C * e
    _stress[_qp](0,0) = _lambda_o * ( _mechanical_strain[_qp](0,0) + _mechanical_strain[_qp](1,1) + _mechanical_strain[_qp](2,2) ) + 2 * shear_modulus * _mechanical_strain[_qp](0,0);
    _stress[_qp](1,1) = _lambda_o * ( _mechanical_strain[_qp](0,0) + _mechanical_strain[_qp](1,1) + _mechanical_strain[_qp](2,2) ) + 2 * shear_modulus * _mechanical_strain[_qp](1,1);
    _stress[_qp](2,2) = _lambda_o * ( _mechanical_strain[_qp](0,0) + _mechanical_strain[_qp](1,1) + _mechanical_strain[_qp](2,2) ) + 2 * shear_modulus * _mechanical_strain[_qp](2,2);
    _stress[_qp](0,1) = 2 * shear_modulus * _mechanical_strain[_qp](0,1); _stress[_qp](1,0) = 2 * shear_modulus * _mechanical_strain[_qp](1,0);
    _stress[_qp](0,2) = 2 * shear_modulus * _mechanical_strain[_qp](0,2); _stress[_qp](2,0) = 2 * shear_modulus * _mechanical_strain[_qp](2,0);
    _stress[_qp](1,2) = 2 * shear_modulus * _mechanical_strain[_qp](1,2); _stress[_qp](2,1) = 2 * shear_modulus * _mechanical_strain[_qp](2,1);

    // Assign value for elastic strain, which is equal to the mechanical strain
    _elastic_strain[_qp] = _mechanical_strain[_qp];

  }
  else{

    /* compute alpha */
    //compute forcing term
    ADReal alpha_forcingterm;
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
    ADReal alpha_damagedvar = _alpha_damagedvar_old[_qp] + _dt * alpha_forcingterm;

    //check alpha within range
    if ( alpha_damagedvar < 0 ){ alpha_damagedvar = 0.0; }
    else if ( alpha_damagedvar > 1 ){ alpha_damagedvar = 1.0; }
    else{} 

    //check below initial damage (fix initial damage)
    if ( alpha_damagedvar < _initial_damage[_qp] ){ alpha_damagedvar = _initial_damage[_qp]; }
    else{}

    _alpha_damagedvar[_qp] = alpha_damagedvar;

    /* compute B */
    ADReal C_B = _CdCb_multiplier * _Cd_constant;

    //alphacr function
    ADReal alphacr;
    if ( _xi_old[_qp] < _xi_0 ){ alphacr = 1.0;} 
    else if ( _xi_old[_qp] > _xi_0 && _xi_old[_qp] <= _xi_1 ){ alphacr = ((_xi_old[_qp]*2.76e5-7.100521107637101e2*_xi_old[_qp]*7.5e2-7.100521107637101e2*1.4e3-7.100521107637101e2*std::pow(_xi_old[_qp],3)*1.25e2+std::pow(_xi_old[_qp],3)*4.6e4+std::sqrt((7.100521107637101e2*3.68e2-3.19799e5)*(_xi_old[_qp]*-1.44e3-std::pow(_xi_old[_qp],2)*2.1e3+std::pow(_xi_old[_qp],3)*5.6e2+std::pow(_xi_old[_qp],4)*3.0e2+std::pow(_xi_old[_qp],6)*2.5e1+3.576e3)*(-3.590922148807814e-1))*5.9e1+5.152e5)*(-5.9e1/4.0))/(_xi_old[_qp]*3.837588e7-7.100521107637101e2*_xi_old[_qp]*4.416e4+7.100521107637101e2*4.048e3-7.100521107637101e2*std::pow(_xi_old[_qp],2)*2.76e4+std::pow(_xi_old[_qp],2)*2.3984925e7-3.517789e6);}
    else if ( _xi_old[_qp] > _xi_1 && _xi_old[_qp] <= _xi_max ){ alphacr = 6.0e10/(7.100521107637101e2*1.627118644067797e8+_xi_old[_qp]*(7.100521107637101e2*1.016949152542373e8-3.742372881355932e10)-5.987796610169492e10);}
    else{std::cout<<"xi: "<<_xi_old[_qp]<<std::endl;mooseError("xi exceeds the maximum allowable range!");}

    //compute forcing func
    ADReal Prob = 1.0 / ( std::exp( (alphacr - _alpha_damagedvar_old[_qp]) / _beta_width ) + 1.0 );
    ADReal B_forcingterm;
    if ( _xi_old[_qp] >= _xi_d && _xi_old[_qp] <= _xi_max ){
      B_forcingterm = 1.0 * C_B * Prob * (1-_B_breakagevar_old[_qp]) * _I2_old[_qp] * (_xi_old[_qp] - _xi_d); //could heal if xi < xi_0
    }
    else if ( _xi_old[_qp] < _xi_d && _xi_old[_qp] >= _xi_min ){
      B_forcingterm = 0.0; //1.0 * _CBH_constant * _I2_old[_qp] * ( _xi_old[_qp] - _xi_d ); //close healing
    }
    else{
      mooseError("xi_old is OUT-OF-RANGE!.");
    }

    ADReal B_damagedvar = _B_breakagevar_old[_qp] + _dt * B_forcingterm;

    //check breakage within range
    if ( B_damagedvar < 0 ){ B_damagedvar = 0.0; }
    else if ( B_damagedvar > 1 ){ B_damagedvar = 1.0; }
    else{}   

    _B_breakagevar[_qp] = B_damagedvar;

    /* update modulus */
    ADReal shear_modulus = _shear_modulus_o + _alpha_damagedvar_old[_qp] * _xi_0 * _gamma_damaged_r;
    ADReal gamma_damaged = _alpha_damagedvar_old[_qp] * _gamma_damaged_r;
    _shear_modulus[_qp] = shear_modulus;
    _gamma_damaged[_qp] = gamma_damaged;

    /* compute strain */
    ADRankTwoTensor eps_p = _eps_p_old[_qp] + _dt * _C_g * std::pow(_B_breakagevar_old[_qp],_m1) * _sigma_d_old[_qp];
    ADRankTwoTensor eps_e = _mechanical_strain[_qp] - eps_p;
    ADReal I1 = eps_e(0,0) + eps_e(1,1) + eps_e(2,2);
    ADReal I2 = eps_e(0,0) * eps_e(0,0) + eps_e(1,1) * eps_e(1,1) + eps_e(2,2) * eps_e(2,2) + 2 * eps_e(0,1) * eps_e(0,1) + 2 * eps_e(0,2) * eps_e(0,2) + 2 * eps_e(1,2) * eps_e(1,2);
    ADReal xi = I1/std::sqrt(I2);

    //Represent sigma (solid(s) + granular(b))
    ADRankTwoTensor sigma_s;
    ADRankTwoTensor sigma_b;
    ADRankTwoTensor sigma_total;
    ADRankTwoTensor sigma_d;
    const auto I = ADRankTwoTensor::Identity();

    sigma_s(0,0) = ( _lambda_o - gamma_damaged / xi ) * I1 + ( 2 * shear_modulus - gamma_damaged * xi ) * eps_e(0,0);
    sigma_b(0,0) = ( 2 * _a2 + _a1 / xi + 3 * _a3 * xi ) * I1 + ( 2 * _a0 + _a1 * xi - _a3 * std::pow(xi,3) ) * eps_e(0,0);
    sigma_s(1,1) = ( _lambda_o - gamma_damaged / xi ) * I1 + ( 2 * shear_modulus - gamma_damaged * xi ) * eps_e(1,1);
    sigma_b(1,1) = ( 2 * _a2 + _a1 / xi + 3 * _a3 * xi ) * I1 + ( 2 * _a0 + _a1 * xi - _a3 * std::pow(xi,3) ) * eps_e(1,1);
    sigma_s(2,2) = ( _lambda_o - gamma_damaged / xi ) * I1 + ( 2 * shear_modulus - gamma_damaged * xi ) * eps_e(2,2);
    sigma_b(2,2) = ( 2 * _a2 + _a1 / xi + 3 * _a3 * xi ) * I1 + ( 2 * _a0 + _a1 * xi - _a3 * std::pow(xi,3) ) * eps_e(2,2);
    sigma_s(0,1)  = ( 2 * shear_modulus - gamma_damaged * xi    ) * eps_e(0,1); sigma_s(1,0)  = ( 2 * shear_modulus - gamma_damaged * xi    ) * eps_e(1,0);
    sigma_b(0,1)  = ( 2 * _a0 + _a1 * xi - _a3 * std::pow(xi,3) ) * eps_e(0,1); sigma_b(1,0)  = ( 2 * _a0 + _a1 * xi - _a3 * std::pow(xi,3) ) * eps_e(1,0);
    sigma_s(0,2)  = ( 2 * shear_modulus - gamma_damaged * xi    ) * eps_e(0,2); sigma_s(2,0)  = ( 2 * _a0 + _a1 * xi - _a3 * std::pow(xi,3) ) * eps_e(2,0);
    sigma_b(0,2)  = ( 2 * _a0 + _a1 * xi - _a3 * std::pow(xi,3) ) * eps_e(0,2); sigma_b(2,0)  = ( 2 * _a0 + _a1 * xi - _a3 * std::pow(xi,3) ) * eps_e(2,0);
    sigma_s(1,2)  = ( 2 * shear_modulus - gamma_damaged * xi    ) * eps_e(1,2); sigma_s(2,1)  = ( 2 * _a0 + _a1 * xi - _a3 * std::pow(xi,3) ) * eps_e(2,1);
    sigma_b(1,2)  = ( 2 * _a0 + _a1 * xi - _a3 * std::pow(xi,3) ) * eps_e(1,2); sigma_b(2,1)  = ( 2 * _a0 + _a1 * xi - _a3 * std::pow(xi,3) ) * eps_e(2,1);
    
    sigma_total = (1 - _B_breakagevar_old[_qp]) * sigma_s + _B_breakagevar_old[_qp] * sigma_b;
    sigma_d = sigma_total - 1/3 * (sigma_total(0,0) + sigma_total(1,1) + sigma_total(2,2)) * I;

    _eps_p[_qp] = eps_p;
    _eps_e[_qp] = eps_e;
    _I1[_qp] = I1;
    _I2[_qp] = I2;
    _xi[_qp] = xi;
    _sigma_d[_qp] = sigma_d;

    // Rotate the stress state to the current configuration
    _stress[_qp] = sigma_total;

    // Assign value for elastic strain, which is equal to the mechanical strain
    _elastic_strain[_qp] = _mechanical_strain[_qp];
  }
}