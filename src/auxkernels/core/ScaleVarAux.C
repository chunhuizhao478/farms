#include "ScaleVarAux.h"

registerMooseObject("farmsApp", ScaleVarAux);

InputParameters
ScaleVarAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addRequiredCoupledVar("coupled","Nonlinear Variable that needed to be scaled");
  params.addRequiredCoupledVar("scale","scale to apply");
  return params;
}

ScaleVarAux::ScaleVarAux(const InputParameters & parameters)
  : AuxKernel(parameters),
  
  _coupled_val(coupledValue("coupled")),
  _scale(coupledValue("scale"))

{
}

Real
ScaleVarAux::computeValue()
{
  return 1.0 * _coupled_val[_qp] / _scale[_qp];
}