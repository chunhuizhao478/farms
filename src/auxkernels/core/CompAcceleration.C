/*
AuxKernel of Passing Variable Time Derivative 
*/

#include "CompAcceleration.h"

registerMooseObject("farmsApp", CompAcceleration);

InputParameters
CompAcceleration::validParams()
{
  InputParameters params = AuxKernel::validParams();

  params.addRequiredCoupledVar("coupled","Nonlinear Variable that needed to be taken time derivative of");

  return params;
}

CompAcceleration::CompAcceleration(const InputParameters & parameters)
  : AuxKernel(parameters),
  
  //Compute the time derivative of the given variable using "coupledDotDot"
  _coupled_val(coupledDotDot("coupled"))

{
}

Real
CompAcceleration::computeValue()
{
  return _coupled_val[_qp];
}