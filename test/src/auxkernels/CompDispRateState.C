#include "CompDispRateState.h"

registerMooseObject("farmsApp", CompDispRateState);

InputParameters
CompDispRateState::validParams()
{
  InputParameters params = AuxKernel::validParams();

  params.addRequiredCoupledVar("current_disp","Current Displacement");
  params.addRequiredCoupledVar("vel_const","Background Constant Velocity");
  return params;
}

CompDispRateState::CompDispRateState(const InputParameters & parameters)
  : AuxKernel(parameters),
  
  _current_disp(coupledValue("current_disp")),
  _vel_const(coupledValue("vel_const"))

{
}

Real
CompDispRateState::computeValue()
{ 
  return _current_disp[_qp] + ( _vel_const[_qp] ) * _t;
}