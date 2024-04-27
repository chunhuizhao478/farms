#include "CheckAlphaBorehole.h"

registerMooseObject("farmsApp", CheckAlphaBorehole);

InputParameters
CheckAlphaBorehole::validParams()
{
  InputParameters params = AuxKernel::validParams();

  params.addRequiredCoupledVar("coupled","Either Alpha or B that needs to be checked");

  return params;
}

CheckAlphaBorehole::CheckAlphaBorehole(const InputParameters & parameters)
  : AuxKernel(parameters),
  
  _coupled_val(coupledValue("coupled"))

{
}

Real
CheckAlphaBorehole::computeValue()
{
  
  Real var = _coupled_val[_qp];
  Real var_out;
  if ( var < 0 )
  {
    var_out = 0;
  }
  else if ( var > 0.6 )
  {
    var_out = 0.6;
  } 
  else{
    var_out = var;
  }
  return var_out;
}