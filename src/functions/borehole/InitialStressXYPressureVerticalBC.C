/*
Define Function for Initial Shear Stress Hetergeneoity in XY Direction 
Will be Passed into AuxVariables
Chunhui Zhao
*/

#include "InitialStressXYPressureVerticalBC.h"

registerMooseObject("farmsApp", InitialStressXYPressureVerticalBC);

InputParameters
InitialStressXYPressureVerticalBC::validParams()
{
  InputParameters params = Function::validParams();
  return params;
}

InitialStressXYPressureVerticalBC::InitialStressXYPressureVerticalBC(const InputParameters & parameters)
  : Function(parameters)
{
}

Real
InitialStressXYPressureVerticalBC::value(Real t, const Point & /*p*/) const
{

  //pressure magnitude as function of time
  Real p_l = 0.0;
  Real rate1 = 20e6 / (200*5.9e-8); //200 time steps
  Real rate2 = 20e6 / (400*5.9e-8); //400 time steps
  
  //assume linearly ramp up
  Real t1 = 200*5.9e-8;
  if ( t < t1 ){
    p_l = rate1 * t;
  }
  else{
    p_l = rate1 * t1 + rate2 * ( t - t1 );
  }

  //set & check threshold 
  Real threshold = 40e6; 
  
  if ( p_l > threshold ){
    p_l = threshold;
  }else{}

  return p_l;

}