#include "FDCompVarRatev2.h"

registerMooseObject("farmsApp", FDCompVarRatev2);

InputParameters
FDCompVarRatev2::validParams()
{
  InputParameters params = AuxKernel::validParams();

  params.addRequiredCoupledVar("coupled","Nonlinear Variable that needed to be taken time derivative of");

  return params;
}

FDCompVarRatev2::FDCompVarRatev2(const InputParameters & parameters)
  : AuxKernel(parameters),
  
  //Compute the time derivative of the given variable using "coupledDot"
  _coupled_val(coupledValue("coupled")),
  //use "old"
  _coupled_val_old(coupledValueOld("coupled"))

{
}

Real
FDCompVarRatev2::computeValue()
{
  return ( _coupled_val[_qp] - _coupled_val_old[_qp] ) / _dt;
}