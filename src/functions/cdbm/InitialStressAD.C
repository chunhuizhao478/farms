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

  if (y_coord >= 0-1*0.01 and y_coord <= 0){
    if (x_coord >= -0.5-2*0.01 and x_coord <= -0.5+2*0.01){
        T1_o = 40e6;
    }
    else if (x_coord <= -0.6 || x_coord >= 0.6){
        T1_o = 20e6;
    }
    else{
        T1_o = 20e6;
    }
  }
  else{
    T1_o = 20e6;
  }
  
  return T1_o;

}