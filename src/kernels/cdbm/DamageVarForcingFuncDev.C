/*
Implementation of Damage Evolution Forcing Function (F) :

Strong Form: 

d alpha / dt = F
F = (1-B)[Cd I_2 (xi - xi_o) + D grad^2 alpha] if xi >= xi_o
F = (1-B)[C1 exp(alpha/C2) I2 (xi - xi_o) + D grad^2 alpha] if xi <= xi_o

Weak Form:

int( d(alpha)/dt * v ) - int( (1-B) (Cd I2 (xi - xi_o) * v + D d(alpha)/dx dv/dx ) = 0 if xi >= xi_o

int( d(alpha)/dt * v ) - int( (1-B) (C1 exp(alpha/C2) I2 (xi - xi_o) * v + D d(alpha)/dx dv/dx ) = 0 if xi <= xi_o

- 10/5/2023 Chunhui Zhao

Include power-law correction on Cd (function of strain rate e)
if (e < 1e-4){ Cd = Cd_min };
else{ Cd = pow(10, log10(1+m*log10(e/1e-4)) ) }

*/

#include "DamageVarForcingFuncDev.h"

registerMooseObject("farmsApp", DamageVarForcingFuncDev);

InputParameters
DamageVarForcingFuncDev::validParams()
{
  InputParameters params = Kernel::validParams();

  //constant parameters
  params.addParam<Real>(   "C_d_min", 10.0, "coefficient gives positive damage evolution (small strain e < 1e-4 threshold value)");
  params.addRequiredParam<Real>(     "D", "coefficient gives diffusion magnitude of damage evolution");
  params.addRequiredParam<Real>(   "C_1", "coefficient of healing for damage evolution");
  params.addRequiredParam<Real>(   "C_2", "coefficient of healing for damage evolution");
  params.addRequiredParam<Real>(  "xi_0", "strain invariants ratio: onset of damage evolution");
  params.addRequiredParam<Real>("xi_min", "strain invariant ratio at minimum value");
  params.addRequiredParam<Real>("xi_max", "strain invariant ratio at maximum value");
  params.addRequiredParam<Real>(     "m", "Cd power-law correction index");
  params.addParam<Real>("mechanical_strain_rate_threshold", 0, "threshold value for strain rate such that Cd takes constant value Cd_min if strain rate below this value.");
  params.addParam<Real>( "scale", 1.0, "scale the Cd power-law");

  //variable parameters
  params.addRequiredCoupledVar("alpha_old", "damage variable at previous time step");
  params.addRequiredCoupledVar(    "B_old", "breakage variable at previous time step");
  params.addRequiredCoupledVar(   "xi_old", "strain invariant ratio at previous time step");
  params.addRequiredCoupledVar(   "I2_old", "second strain invariant at previous time step");
  params.addCoupledVar("mechanical_strain_rate", 0.0, "strain rate");

  //add options
  params.addRequiredParam<int>( "option", "option 1 : Cd power-law; option 2 : use constant Cd");
  params.addParam<Real>( "Cd_constant", 0.0, "constant Cd value for option 2 only");

  //fix diffuison
  params.addRequiredParam<Real>("shear_modulus_o", "initial shear modulus");
  params.addRequiredParam<Real>("lambda_o", "initial lame constant");

  //add healing
  params.addParam<bool>("healing", false, "if turn on healing, true = on, false = off, default is false = off");

  return params;
}

DamageVarForcingFuncDev::DamageVarForcingFuncDev(const InputParameters & parameters)
 : Kernel(parameters),
  _Cd_min(getParam<Real>("C_d_min")),
  _D(getParam<Real>("D")),
  _C1(getParam<Real>("C_1")),
  _C2(getParam<Real>("C_2")),
  _xi_0(getParam<Real>("xi_0")),
  _xi_min(getParam<Real>("xi_min")),
  _xi_max(getParam<Real>("xi_max")),
  _scale(getParam<Real>("scale")),
  _m(getParam<Real>("m")),
  _mechanical_strain_rate_threshold(getParam<Real>("mechanical_strain_rate_threshold")),
  _alpha_old(coupledValue("alpha_old")),
  _B_old(coupledValue("B_old")),
  _xi_old(coupledValue("xi_old")),
  _I2_old(coupledValue("I2_old")),
  _mechanical_strain_rate(coupledValue("mechanical_strain_rate")),
  _option(getParam<int>("option")),
  _Cd_constant(getParam<Real>("Cd_constant")),
  _shear_modulus_o(getParam<Real>("shear_modulus_o")),
  _lambda_o(getParam<Real>("lambda_o")),
  _healing(getParam<bool>("healing"))
{
}

Real
DamageVarForcingFuncDev::computeQpResidual()
{ 
  
  //Power-law correction
  //Initialize Cd
  Real Cd = 0;
  //Check options
  if ( _option == 1 ){

    //close _option == 1
    mooseError("Option 1 is NOT available!");

    //power-law correction on coefficient Cd(function of strain rate)
    if ( _mechanical_strain_rate[_qp] < _mechanical_strain_rate_threshold ) //Cd follows power-law
    {
      Cd = _scale * pow(10, 1 + _m * log10( _mechanical_strain_rate[_qp] / _mechanical_strain_rate_threshold ) ) * _Cd_min;
    }
    else if ( _mechanical_strain_rate[_qp] < 0 ){ //Cd remains constant
      Cd = _Cd_min;
    }

  }
  else if ( _option == 2 ){

    if ( _Cd_constant == 0.0 ){
      mooseError("For option 2, need to provide nonzero Cd_constant value !");
    }
    else{
      Cd = _Cd_constant;
    }

  }
  else{
    mooseError("Please provide valid option number!");
  }

  //Compute Diffusion coefficient
  Real YoungE = _shear_modulus_o * ( 3 * _lambda_o + 2 * _shear_modulus_o ) / ( _lambda_o + _shear_modulus_o );
  Real Diffusion_Coeff = _D * Cd / YoungE;

  //weak form for damage variable evolution
  if ( _xi_old[_qp] >= _xi_0 && _xi_old[_qp] <= _xi_max ){
    return -1 * (1 - _B_old[_qp]) * ( Cd * _I2_old[_qp] * ( _xi_old[_qp] - _xi_0 ) * _test[_i][_qp] + Diffusion_Coeff * _grad_u[_qp] * _grad_test[_i][_qp] );
  }
  else if ( _xi_old[_qp] < _xi_0 && _xi_old[_qp] >= _xi_min ){
    //with healing
    if ( _healing == true ){
      return -1 * (1 - _B_old[_qp]) * ( _C1 * exp(_alpha_old[_qp]/_C2) * _I2_old[_qp] * ( _xi_old[_qp] - _xi_0 ) * _test[_i][_qp] + Diffusion_Coeff * _grad_u[_qp] * _grad_test[_i][_qp] );
    }
    //no healing
    else{
      return 0.0;
    }
  }
  else{
    mooseError("xi_old is OUT-OF-RANGE!.");
    return 0;
  }
}

Real
DamageVarForcingFuncDev::computeQpJacobian()
{
  return 0.0;
}