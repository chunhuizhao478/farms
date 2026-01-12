#include "FarmsScaleVarAux.h"
#include "Assembly.h"

registerMooseObject("farmsApp", FarmsScaleVarAux);

InputParameters
FarmsScaleVarAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addRequiredCoupledVar("coupled","Nonlinear Variable that needed to be scaled");
  return params;
}

FarmsScaleVarAux::FarmsScaleVarAux(const InputParameters & parameters)
  : AuxKernel(parameters),
  
  _coupled_val(coupledValue("coupled")),
  _current_elem_volume(_assembly.elemVolume()),
  _current_side_volume(_assembly.sideElemVolume())

{
}

Real
FarmsScaleVarAux::computeValue()
{
  return 1.0 * _coupled_val[_qp] / 40000;
}