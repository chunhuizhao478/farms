#include "CheckAlphaB.h"

registerMooseObject("farmsApp", CheckAlphaB);

InputParameters
CheckAlphaB::validParams()
{
  InputParameters params = AuxKernel::validParams();

  params.addRequiredCoupledVar("coupled","Either Alpha or B that needs to be checked");

  return params;
}

CheckAlphaB::CheckAlphaB(const InputParameters & parameters)
  : AuxKernel(parameters),
  
  _coupled_val(coupledValue("coupled"))

{
}

Real
CheckAlphaB::computeValue()
{
  
  Real var = _coupled_val[_qp];
  Real var_out;
  if ( var < 0 )
  {
    var_out = 0;
  }
  else if ( var > 1 )
  {
    var_out = 1;
  } 
  else{
    var_out = var;
  }
  return var_out;
}