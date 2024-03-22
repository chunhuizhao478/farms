/*
Implementation of Damage Evolution Forcing Function (F) :

Strong Form: 

d alpha / dt = F
F = (1-B)[Cd I_2 (xi - xi_o) + D grad^2 alpha] if xi >= xi_o
F = (1-B)[C1 exp(alpha/C2) I2 (xi - xi_o) + D grad^2 alpha] if xi <= xi_o

Weak Form:

int( d(alpha)/dt * v ) - int( (1-B) (Cd I2 (xi - xi_o) * v + D d(alpha)/dx dv/dx ) = 0 if xi >= xi_o

int( d(alpha)/dt * v ) - int( (1-B) (C1 exp(alpha/C2) I2 (xi - xi_o) * v + D d(alpha)/dx dv/dx ) = 0 if xi <= xi_o

*/

#include "ADDamageVarForcingFunc.h"

registerMooseObject("farmsApp", ADDamageVarForcingFunc);

InputParameters
ADDamageVarForcingFunc::validParams()
{
  InputParameters params = ADKernel::validParams();

  //constant parameters
  params.addParam<Real>(   "Cd_constant", 1.0, "coefficient gives positive damage evolution");
  params.addParam<Real>(     "D", 1.0, "coefficient gives diffusion magnitude of damage evolution");
  params.addParam<Real>(   "C_1", 1.0, "coefficient of healing for damage evolution");
  params.addParam<Real>(   "C_2", 1.0, "coefficient of healing for damage evolution");
  params.addParam<Real>(  "xi_0", 1.0, "strain invariants ratio: onset of damage evolution");
  params.addParam<Real>("xi_min", 1.0, "strain invariant ratio at minimum value");
  params.addParam<Real>("xi_max", 1.0, "strain invariant ratio at maximum value");

  //variable parameters
  params.addRequiredCoupledVar("alpha_old", "damage variable at previous time step");
  params.addRequiredCoupledVar(    "B_old", "breakage variable at previous time step");
  params.addRequiredCoupledVar(   "xi_old", "strain invariant ratio at previous time step");
  params.addRequiredCoupledVar(   "I2_old", "second strain invariant at previous time step");

  return params;
}

ADDamageVarForcingFunc::ADDamageVarForcingFunc(const InputParameters & parameters)
 : ADKernel(parameters),
  _Cd(getParam<Real>("Cd_constant")),
  _D(getParam<Real>("D")),
  _C1(getParam<Real>("C_1")),
  _C2(getParam<Real>("C_2")),
  _xi_0(getParam<Real>("xi_0")),
  _xi_min(getParam<Real>("xi_min")),
  _xi_max(getParam<Real>("xi_max")),
  _alpha_old(adCoupledValue("alpha_old")),
  _B_old(adCoupledValue("B_old")),
  _xi_old(adCoupledValue("xi_old")),
  _I2_old(adCoupledValue("I2_old"))
{
}

ADReal
ADDamageVarForcingFunc::computeQpResidual()
{ 
  // if ( _xi_old[_qp] >= _xi_0 && _xi_old[_qp] <= _xi_max ){
  //   return -1 * (1 - _B_old[_qp]) * ( _Cd * _I2_old[_qp] * ( _xi_old[_qp] - _xi_0 ) * _test[_i][_qp] + _D * _grad_u[_qp] * _grad_test[_i][_qp] );
  // }
  // else if ( _xi_old[_qp] < _xi_0 && _xi_old[_qp] >= _xi_min ){
  //   //with healing
  //   return -1 * (1 - _B_old[_qp]) * ( _C1 * std::exp(_alpha_old[_qp]/_C2) * _I2_old[_qp] * ( _xi_old[_qp] - _xi_0 ) * _test[_i][_qp] + _D * _grad_u[_qp] * _grad_test[_i][_qp] );
  //   //no healing
  //   //return 0.0;
  // }
  // else{
  //   mooseError("xi_old is OUT-OF-RANGE!.");
  //   return 0;
  // }

  if ( _xi_old[_qp] >= _xi_0 && _xi_old[_qp] <= _xi_max ){
    return -1 * (1 - _B_old[_qp]) * ( _Cd * _I2_old[_qp] * ( _xi_old[_qp] - _xi_0 ) );
  }
  else if ( _xi_old[_qp] < _xi_0 && _xi_old[_qp] >= _xi_min ){
    //with healing
    return -1 * (1 - _B_old[_qp]) * ( _C1 * std::exp(_alpha_old[_qp]/_C2) * _I2_old[_qp] * ( _xi_old[_qp] - _xi_0 ) );
    //no healing
    //return 0.0;
  }
  else{
    mooseError("xi_old is OUT-OF-RANGE!.");
    return 0;
  }

  //ggw183
  // if ( _xi_old[_qp] >= _xi_0 && _xi_old[_qp] <= _xi_max ){
  //   return -1 * ( _Cd * _I2_old[_qp] * ( _xi_old[_qp] - _xi_0 ) * _test[_i][_qp] + _D * _grad_u[_qp] * _grad_test[_i][_qp] );
  // }
  // else if ( _xi_old[_qp] < _xi_0 && _xi_old[_qp] >= _xi_min ){
  //   //with healing
  //   return -1 * ( _C1 * std::exp(_alpha_old[_qp]/_C2) * _I2_old[_qp] * ( _xi_old[_qp] - _xi_0 ) * _test[_i][_qp] + _D * _grad_u[_qp] * _grad_test[_i][_qp] );
  //   //no healing
  //   //return 0.0;
  // }
  // else{
  //   mooseError("xi_old is OUT-OF-RANGE!.");
  //   return 0;
  // }
}
