#include "InitialStressAD.h"

registerMooseObject("farmsApp", InitialStressAD);

InputParameters
InitialStressAD::validParams()
{
  InputParameters params = Function::validParams();
  return params;
}

InitialStressAD::InitialStressAD(const InputParameters & parameters)
  : Function(parameters)
{
}

Real
InitialStressAD::value(Real /*t*/, const Point & p) const
{


  Real x_coord = p(0); //along the strike direction
  Real y_coord = p(1); //along the normal direction

  Real T1_o = 0;

  if (y_coord >= 0-2*0.01 and y_coord <= 0+2*0.01){
    if (x_coord >= -1.0-2*0.01 and x_coord <= -1.0+2*0.01){
        T1_o = 35e6;
    }
    else if (x_coord <= -1.1 || x_coord >= 1.1){
        T1_o = 30e6;
    }
    else{
        T1_o = 30e6;
    }
  }
  else{
    T1_o = 30e6;
  }
  
  return T1_o;

}