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

  //pressure magnitude as function of time
  Real p_l = 0;
  Real rate = 15e6 / (200*5.9e-8);
  
  //assume linearly ramp up
  p_l = rate * t;

  //set & check threshold 
  Real threshold = 15e6;
  
  if ( p_l > threshold ){
    p_l = threshold;
  }else{}

  return p_l;

}