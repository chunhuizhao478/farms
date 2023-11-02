/*
Define Function for Initial Shear Stress Hetergeneoity in XY Direction 
Will be Passed into AuxVariables
Chunhui Zhao
*/

#include "InitialStressXYnetwork_old.h"

registerMooseObject("farmsApp", InitialStressXYnetwork_old);

InputParameters
InitialStressXYnetwork_old::validParams()
{
  InputParameters params = Function::validParams();
  return params;
}

InitialStressXYnetwork_old::InitialStressXYnetwork_old(const InputParameters & parameters)
  : Function(parameters)
{
}

Real
InitialStressXYnetwork_old::value(Real /*t*/, const Point & p) const
{

  Real x_coord = p(0); //along the strike direction
  Real y_coord = p(1); //along the normal direction

  Real T1_o = 0;

  Real origin_x = 208.7;
  Real origin_y = -1345.82;

  if (x_coord >= origin_x-300 and x_coord <= origin_x+300 and y_coord >= origin_y-300 and y_coord <= origin_y+300)
  {
    T1_o = 88.2e6; //1% overstress
    // T1_o = 96e6;
  }
  else{
    T1_o = 70e6;
  }
  
  // return T1_o;

  // Real T1_o = 0;

  // if (x_coord >= -400 and x_coord <= 400 and y_coord >= -400 and y_coord <= 400)
  // {
  //   T1_o = 81.6e6;
  // }
  // else{
  //   T1_o = 70e6;
  // }
  
  return T1_o;

}