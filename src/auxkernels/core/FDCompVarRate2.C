#include "FDCompVarRate2.h"

registerMooseObject("farmsApp", FDCompVarRate2);

InputParameters
FDCompVarRate2::validParams()
{
  InputParameters params = AuxKernel::validParams();

  params.addRequiredCoupledVar("coupled","Nonlinear Variable that needed to be taken time derivative of");

  return params;
}

FDCompVarRate2::FDCompVarRate2(const InputParameters & parameters)
  : AuxKernel(parameters),
  
  //Compute the time derivative of the given variable using "coupledDot"
  _coupled_val(coupledValue("coupled")),
  _coupled_val_old(coupledValueOld("coupled"))

{
}

Real
FDCompVarRate2::computeValue()
{
  return ( _coupled_val[_qp] - _coupled_val_old[_qp] ) / _dt;
}