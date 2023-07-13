#include "CompTJumpRateState.h"

registerMooseObject("farmsApp", CompTJumpRateState);

InputParameters
CompTJumpRateState::validParams()
{
  InputParameters params = AuxKernel::validParams();

  params.addRequiredCoupledVar("coupled","Nonlinear Variable that needed to be taken time derivative of");
  params.addParam<Real>("sliprate_bd",0,"background slip rate");
  return params;
}

CompTJumpRateState::CompTJumpRateState(const InputParameters & parameters)
  : AuxKernel(parameters),
  
  _coupled_val(coupledValue("coupled")),
  _sliprate_bd(getParam<Real>("sliprate_bd"))

{
}

Real
CompTJumpRateState::computeValue()
{
  return _coupled_val[_qp] + _sliprate_bd * _dt;
}