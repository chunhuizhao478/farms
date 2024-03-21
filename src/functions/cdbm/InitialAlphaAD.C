#include "InitialAlphaAD.h"

registerMooseObject("farmsApp", InitialAlphaAD);

InputParameters
InitialAlphaAD::validParams()
{
  InputParameters params = Function::validParams();
  return params;
}

InitialAlphaAD::InitialAlphaAD(const InputParameters & parameters)
  : Function(parameters)
{
}

Real
InitialAlphaAD::value(Real /*t*/, const Point & p) const
{


  Real x_coord = p(0); //along the strike direction
  Real y_coord = p(1); //along the normal direction

  Real alpha_o = 0;

  if (y_coord >= 0-2*0.01 and y_coord <= 0+2*0.01){
    if (x_coord >= -1.0-2*0.01 and x_coord <= -1.0+2*0.01){
        alpha_o = 0.8;
    }
    else{
        alpha_o = 0.7;
    }
  }
  else{
    alpha_o = 0.0;
  }
  
  return alpha_o;

}