#include "CheckAlphaCycleSim.h"

registerMooseObject("farmsApp", CheckAlphaCycleSim);

InputParameters
CheckAlphaCycleSim::validParams()
{
  InputParameters params = AuxKernel::validParams();

  params.addRequiredCoupledVar("coupled","Either Alpha or B that needs to be checked");

  return params;
}

CheckAlphaCycleSim::CheckAlphaCycleSim(const InputParameters & parameters)
  : AuxKernel(parameters),
  
  _coupled_val(coupledValue("coupled"))

{
}

Real
CheckAlphaCycleSim::computeValue()
{

  //hardcode initial condition
  Real x_coord = _q_point[_qp](0); //along the strike direction
  Real y_coord = _q_point[_qp](1); //along the normal direction
  Real initial_damage = 0.0;

  if (y_coord >= 0-1*0.1 and y_coord <= 0+1*0.1){
    if (x_coord >= 0-1*0.1 and x_coord <= 0+1*0.1){
        initial_damage = 0.7;
    }
    else if(x_coord <= -5 || x_coord >= 5){
        initial_damage = 0.7;
    }
    else{
        initial_damage = 0.7;
    }
  }
  else{
    initial_damage = 0.0;
  }
  //hardcode initial condition

  Real var = _coupled_val[_qp]; 
  Real var_out;

  if ( var < initial_damage ){
    return initial_damage;
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