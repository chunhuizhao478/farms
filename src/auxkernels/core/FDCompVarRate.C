#include "FDCompVarRate.h"

registerMooseObject("farmsApp", FDCompVarRate);

InputParameters
FDCompVarRate::validParams()
{
  InputParameters params = AuxKernel::validParams();

  params.addRequiredCoupledVar("coupled","Nonlinear Variable that needed to be taken time derivative of");

  return params;
}

FDCompVarRate::FDCompVarRate(const InputParameters & parameters)
  : AuxKernel(parameters),
  
  //Compute the time derivative of the given variable using "coupledDot"
  _coupled_val(coupledValue("coupled")),
  _coupled_val_old(coupledValueOlder("coupled"))

{
}

Real
FDCompVarRate::computeValue()
{
  return ( _coupled_val[_qp] - _coupled_val_old[_qp] ) / _dt;
}