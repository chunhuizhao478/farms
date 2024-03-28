#include "CheckAlphaTDZ.h"

registerMooseObject("farmsApp", CheckAlphaTDZ);

InputParameters
CheckAlphaTDZ::validParams()
{
  InputParameters params = AuxKernel::validParams();

  params.addRequiredCoupledVar("coupled","Either Alpha or B that needs to be checked");
  params.addRequiredCoupledVar("initial_damage","Initial damage distribution");

  return params;
}

CheckAlphaTDZ::CheckAlphaTDZ(const InputParameters & parameters)
  : AuxKernel(parameters),
  
  _coupled_val(coupledValue("coupled")),

  _initial_damage(coupledValue("initial_damage"))

{
}

Real
CheckAlphaTDZ::computeValue()
{
  
  Real var = _coupled_val[_qp]; 
  Real var_out;

  if ( var < _initial_damage[_qp] ){
    return _initial_damage[_qp];
  }
  else{
    if ( var < 0 )
    {
      var_out = 0;
    }
    else if ( var > 1.0 )
    {
      var_out = 1.0;
    } 
    else{
      var_out = var;
    } 
    return var_out;
  }
}