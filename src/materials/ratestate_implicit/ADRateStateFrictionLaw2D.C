#include "ADRateStateFrictionLaw2D.h"

registerMooseObject("farmsApp", ADRateStateFrictionLaw2D);

InputParameters
ADRateStateFrictionLaw2D::validParams()
{
  InputParameters params = ADRateStateFrictionLaw2D::validParams();
  return params;
}

ADRateStateFrictionLaw2D::ADRateStateFrictionLaw2D(const InputParameters & parameters)
  : ADCZMComputeLocalTractionTotalBaseRSF2D(parameters),
  _f_o(getParam<Real>("f_o")),
  _rsf_a(getParam<Real>("rsf_a")),
  _rsf_b(getParam<Real>("rsf_b")),
  _rsf_L(getParam<Real>("rsf_L")),
  _delta_o(getParam<Real>("delta_o")),
  _Vini(getParam<Real>("Vini")),
  _Tn_o(getParam<Real>("Tn_o")),
  _Ts_o(getParam<Real>("Ts_o")),
  _interface_displacement_jump_old(getMaterialPropertyOldByName<RealVectorValue>(_base_name + "interface_displacement_jump")),
  _interface_displacement_jump_older(getMaterialPropertyOlderByName<RealVectorValue>(_base_name + "interface_displacement_jump")),
  _statevar_old(getMaterialPropertyOldByName<Real>("statevar"))
{
}

void
ADRateStateFrictionLaw2D::computeInterfaceTraction()
{
  // //Get initial slip rate value
  // ADRealVectorValue Vini_vec(0.0,_Vini,0.0);

  // //Get normal traction
  // //if (_interface_displacement_jump[_qp](0) < 0){}else{}
  // ADReal Tn = std::abs(_Tn_o);

  // //Approximate slip rate
  // ADRealVectorValue interface_displacement_jump_rate = (1.0/_dt) * (_interface_displacement_jump[_qp] - _interface_displacement_jump_old[_qp]) + Vini_vec;
  // ADRealVectorValue interface_displacement_jump_rate_old = (1.0/_dt) * (_interface_displacement_jump_old[_qp] - _interface_displacement_jump_older[_qp]) + Vini_vec;

  // //Get strike direction at t
  // ADReal sliprate_strike = interface_displacement_jump_rate(1);
  // ADReal sliprate_strike_old = interface_displacement_jump_rate_old(1);

  // //Define state variable at t+dt
  // ADReal statevar = ( _rsf_L / sliprate_strike_old) + ( _statevar_old[_qp] - _rsf_L / sliprate_strike_old ) * std::exp( -1 * sliprate_strike_old * _dt / _rsf_L );
  
  // //Define friction mu at t+dt
  // //!arcsinh(z)= ln(z+sqrt(z^2+1))
  // ADReal mu = _rsf_a * std::asinh( sliprate_strike/(2*_delta_o) * std::exp((_f_o + _rsf_b * std::log(_delta_o * statevar/_rsf_L))/_rsf_a) );
  
  // //Compute shear traction at t+dt
  // ADReal Ts = Tn * mu;

  // //Save results
  // _statevar[_qp] = statevar;
  // _sliprate[_qp] = sliprate_strike;

  // //Feed traction into the system
  // //Assign back traction in CZM
  // ADRealVectorValue traction;

  // // traction(0) = 0.0; 
  // // traction(1) = -Ts+_Ts_o; //"-"
  // // traction(2) = 0.0;

  // _interface_traction[_qp] = traction;
}