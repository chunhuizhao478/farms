/*
Define Function for Initial Shear Stress Hetergeneoity in XY Direction 
Will be Passed into AuxVariables
Chunhui Zhao
*/

#include "InitialStressXYPressureBorehole_x_fast.h"

registerMooseObject("farmsApp", InitialStressXYPressureBorehole_x_fast);

InputParameters
InitialStressXYPressureBorehole_x_fast::validParams()
{
  InputParameters params = Function::validParams();
  return params;
}

InitialStressXYPressureBorehole_x_fast::InitialStressXYPressureBorehole_x_fast(const InputParameters & parameters)
  : Function(parameters)
{
}

Real
InitialStressXYPressureBorehole_x_fast::value(Real t, const Point & p) const
{
  
  //get coordinate of current point
  Real x_coord = p(0);
  Real y_coord = p(1);

  //compute ratio
  Real dx = x_coord - 0.0;
  Real dy = y_coord - 0.0;
  Real dl = std::sqrt(dx*dx+dy*dy);

  //pressure magnitude as function of time
  Real p_l = 0.0;
  Real rate = 20e6;
  
  //assume linearly ramp up
  p_l = rate * t;

  //set & check threshold 
  Real threshold = 20e6; //10% of mean stress maximum

  if ( p_l > threshold ){
    p_l = threshold;
  }else{}

  //get x y magnitude components
  Real p_x = p_l * dx / dl;
  Real p_y = p_l * dy / dl;

  //determine sign based on coordinates
  if (x_coord > 0 && y_coord > 0){
      p_x = p_x; p_y = p_y;
  }
  else if(x_coord > 0 && y_coord < 0){
      p_x = p_x; p_y = -1 * p_y;
  }
  else if(x_coord < 0 && y_coord < 0){
      p_x = -1 * p_x; p_y = -1 * p_y;
  }
  else if(x_coord < 0 && y_coord > 0){
      p_x = -1 * p_x; p_y = p_y;
  }
  else if(x_coord == 0.0 && y_coord > 0){
      p_x = 0; p_y = p_l;
  }
  else if(x_coord == 0.0 && y_coord < 0){
      p_x = 0; p_y = -1 * p_l;
  }
  else if(x_coord > 0.0 && y_coord == 0.0){
      p_x = p_l  ; p_y = 0;
  }
  else if(x_coord < 0.0 && y_coord == 0.0){
      p_x = -1 * p_l ; p_y = 0;
  }

  return p_x;

}