#include "InitialStress3DAD.h"

registerMooseObject("farmsApp", InitialStress3DAD);

InputParameters
InitialStress3DAD::validParams()
{
  InputParameters params = Function::validParams();
  return params;
}

InitialStress3DAD::InitialStress3DAD(const InputParameters & parameters)
  : Function(parameters)
{
}

Real
InitialStress3DAD::value(Real /*t*/, const Point & p) const
{


  Real x_coord = p(0); //along the strike direction
  Real y_coord = p(1); //along the normal direction
  Real z_coord = p(2); //along the dip direction

  Real T1_o = 0;

  if (z_coord >= 0-0.01 and z_coord <= 0){
    if (y_coord >= -0.15 and y_coord <= 0.15){
      if (x_coord >= -0.6-2*0.01 and x_coord <= -0.6+2*0.01 and y_coord >= 0-2*0.01 and y_coord <= 0+2*0.01){
          T1_o = 40e6;
      }
      else if (x_coord <= -0.75 || x_coord >= 0.75){
          T1_o = 25e6;
      }
      else{
          T1_o = 25e6;
      }
    }
    else{
      T1_o = 25e6;
    }
  }
  else{
    T1_o = 25e6;
  }
  
  return T1_o;

}