#include "CheckAlpha.h"

registerMooseObject("farmsApp", CheckAlpha);

InputParameters
CheckAlpha::validParams()
{
  InputParameters params = AuxKernel::validParams();

  params.addRequiredCoupledVar("coupled","Either Alpha or B that needs to be checked");

  return params;
}

CheckAlpha::CheckAlpha(const InputParameters & parameters)
  : AuxKernel(parameters),
  
  _coupled_val(coupledValue("coupled"))

{
}

Real
CheckAlpha::computeValue()
{
  
  Real var = _coupled_val[_qp];
  Real var_out;
  if ( var < 0 )
  {
    var_out = 0;
  }
  else if ( var > 0.75 )
  {
    var_out = 0.75;
  } 
  else{
    var_out = var;
  }
  return var_out;
}