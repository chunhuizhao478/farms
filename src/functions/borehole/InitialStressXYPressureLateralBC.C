/*
Define Function for Initial Shear Stress Hetergeneoity in XY Direction 
Will be Passed into AuxVariables
Chunhui Zhao
*/

#include "InitialStressXYPressureLateralBC.h"

registerMooseObject("farmsApp", InitialStressXYPressureLateralBC);

InputParameters
InitialStressXYPressureLateralBC::validParams()
{
  InputParameters params = Function::validParams();
  return params;
}

InitialStressXYPressureLateralBC::InitialStressXYPressureLateralBC(const InputParameters & parameters)
  : Function(parameters)
{
}

Real
InitialStressXYPressureLateralBC::value(Real t, const Point & p) const
{
  
  //get coordinate of current point
  Real x_coord = p(0);
  Real y_coord = p(1);

  //compute ratio
  Real dx = x_coord - 0.0;
  Real dy = y_coord - 0.0;
  Real dl = std::sqrt(dx*dx+dy*dy);

  //pressure magnitude as function of time
  Real p_l = 16.2e6;
  Real rate = 16.2e6;
  
  //assume linearly ramp up
  p_l = rate * t;

  //set & check threshold 
  Real threshold = 16.2e6;
  
  if ( p_l > threshold ){
    p_l = threshold;
  }else{}

  return p_l;

}