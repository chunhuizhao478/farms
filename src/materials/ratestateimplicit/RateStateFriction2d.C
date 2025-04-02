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
  params.addRequiredParam<Real>("state_variable_ref", "Reference state variable");
  params.addRequiredParam<Real>("length_scale_ref", "Reference length scale");
  params.addRequiredParam<Real>("T1_o", "Background shear traction");
  params.addRequiredParam<Real>("T2_o", "Background normal traction");
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
    _state_variable_ref(getParam<Real>("state_variable_ref")),
    _length_scale_ref(getParam<Real>("length_scale_ref")),
    _T1_o(getParam<Real>("T1_o")),
    _T2_o(getParam<Real>("T2_o"))
{
}

void
RateStateFriction2d::computeInterfaceTractionAndDerivatives()
{

  // Compute slip rate
  _slip_rate[_qp] = (_interface_displacement_jump[_qp] - _interface_displacement_jump_old[_qp]) / _dt;

  // Compute slip rate magnitude
  _slip_rate_magnitude[_qp] = _slip_rate[_qp].norm();

  // Compute state variable
  _state_variable[_qp] = _state_variable_old[_qp] + _dt * ( 1 -  _slip_rate_magnitude[_qp] / _length_scale_ref * _state_variable_old[_qp] );
  
  // Compute interface traction (relative to the background traction)
  _interface_traction[_qp](0) = 0.0;
  _interface_traction[_qp](1) = _a * _T2_o * std::asinh( (_slip_rate_magnitude[_qp] / 2 * _slip_rate_ref) * std::exp((_fo + _b * std::log(_slip_rate_ref * _state_variable[_qp] / _length_scale_ref))/(_a)) ) + T1_o;
  _interface_traction[_qp](2) = 0.0;

  // Compute interface traction derivatives
  _dinterface_traction_djump[_qp] = 0.0;
  
}