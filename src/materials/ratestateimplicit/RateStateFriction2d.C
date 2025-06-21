//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "RateStateFriction2d.h"

registerMooseObject("farmsApp", RateStateFriction2d);

InputParameters
RateStateFriction2d::validParams()
{
  InputParameters params = CZMComputeLocalTractionTotalBase::validParams();
  params.addClassDescription("Pure elastic traction separation law.");
  params.addRequiredParam<Real>("fo", "Friction coefficient");
  params.addRequiredParam<Real>("a", "Rate effect parameter");
  params.addRequiredParam<Real>("b", "State effect parameter");
  params.addRequiredParam<Real>("slip_rate_ref", "Reference slip rate");
  params.addRequiredParam<Real>("state_variable_ini", "initial state variable");
  params.addRequiredParam<Real>("length_scale_ref", "Reference length scale");
  params.addRequiredParam<Real>("T1_o", "Background shear traction");
  params.addRequiredParam<Real>("T2_o", "Background normal traction");
  params.addRequiredParam<Real>("slip_rate_ini", "Initial slip rate");
  params.addRequiredCoupledVar("T1_perturb","shear stress perturbation in strike dir");
  return params;
}

RateStateFriction2d::RateStateFriction2d(const InputParameters & parameters)
  : CZMComputeLocalTractionTotalBase(parameters),
    _state_variable(declareProperty<Real>("state_variable")),
    _state_variable_old(getMaterialPropertyOldByName<Real>("state_variable")),
    _slip_rate(declareProperty<RealVectorValue>("slip_rate")),
    _slip_rate_old(getMaterialPropertyOldByName<RealVectorValue>("slip_rate")),
    _slip_rate_magnitude(declareProperty<Real>("slip_rate_magnitude")),
    _slip_rate_magnitude_old(getMaterialPropertyOld<Real>("slip_rate_magnitude")),
    _interface_displacement_jump_old(
      getMaterialPropertyOld<RealVectorValue>(_base_name + "interface_displacement_jump")),
    _fo(getParam<Real>("fo")),
    _a(getParam<Real>("a")),
    _b(getParam<Real>("b")),
    _slip_rate_ref(getParam<Real>("slip_rate_ref")),
    _state_variable_ini(getParam<Real>("state_variable_ini")),
    _length_scale_ref(getParam<Real>("length_scale_ref")),
    _T1_o(getParam<Real>("T1_o")),
    _T2_o(getParam<Real>("T2_o")),
    _slip_rate_ini(getParam<Real>("slip_rate_ini")),
    _T1_perturb(coupledValue("T1_perturb"))
{
}

void
RateStateFriction2d::initQpStatefulProperties()
{
  _interface_traction[_qp] = 0;
  _state_variable[_qp] = _state_variable_ini;
  _slip_rate[_qp] = 0;
  _slip_rate_magnitude[_qp] = _slip_rate_ini;
}

void
RateStateFriction2d::computeInterfaceTractionAndDerivatives()
{

  // Compute slip rate
  _slip_rate[_qp] = (_interface_displacement_jump[_qp] - _interface_displacement_jump_old[_qp]) / _dt;
  _slip_rate[_qp](1) += _slip_rate_ini;

  // Compute slip rate magnitude
  _slip_rate_magnitude[_qp] = std::sqrt( _slip_rate[_qp](1) * _slip_rate[_qp](1) );

  // Compute state variable
  _state_variable[_qp] = _length_scale_ref * ( _dt + _state_variable_old[_qp] ) / (_slip_rate_magnitude[_qp] * _dt + _length_scale_ref);

  // Compute F(V) and derivatives
  Real F_V = (_fo + _b * std::log(_slip_rate_ref * _state_variable[_qp] / _length_scale_ref))/(_a);

  Real dfvdv = - (_b )/(_a * _state_variable[_qp]) * ( _dt * _length_scale_ref * (_dt + _state_variable_old[_qp]) ) / (_slip_rate_magnitude[_qp] * _dt + _length_scale_ref) / (_slip_rate_magnitude[_qp] * _dt + _length_scale_ref);

  // Compute X
  Real X = _slip_rate_magnitude[_qp] / ( 2 * _slip_rate_ref) * std::exp(F_V);

  // Compute tau
  Real tau = _a * (-_T2_o) * std::asinh(X); // + _T1_perturb[_qp];

  // Compute tau derivative: dtau/dx
  Real dtaudx = _a * (-_T2_o) / std::sqrt( 1 + X * X );

  // Compute tau derivative: dXdV
  Real dxdv = std::exp(F_V) / ( 2 * _slip_rate_ref) * ( 1 + _slip_rate_magnitude[_qp] * dfvdv);

  // Compute tau derivatiave: dV/djump
  Real dvdjump = (_slip_rate[_qp](1))/(_slip_rate_magnitude[_qp] * _dt);

  //std::cout << -tau <<std::endl;

  // Compute interface traction (relative to the background traction)
  _interface_traction[_qp](0) = 0.0;
  _interface_traction[_qp](1) = -tau;
  _interface_traction[_qp](2) = 0.0;

  // Compute interface traction derivatives
  _dinterface_traction_djump[_qp] = 0.0;
  _dinterface_traction_djump[_qp](1,1) = dtaudx * dxdv * dvdjump;
  
}