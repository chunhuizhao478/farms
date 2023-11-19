#include "InitialStressXYfrankcase.h"

registerMooseObject("farmsApp", InitialStressXYfrankcase);

InputParameters
InitialStressXYfrankcase::validParams()
{
  InputParameters params = Function::validParams();
  return params;
}

InitialStressXYfrankcase::InitialStressXYfrankcase(const InputParameters & parameters)
  : Function(parameters)
{
}

Real
InitialStressXYfrankcase::value(Real /*t*/, const Point & p) const
{

  Real x_coord = p(0); //along the strike direction
  Real y_coord = p(1); //along the normal direction

  Real T1_o = 0;

  Real center_left_x = -400; Real center_left_y = 40;
  Real center_right_x = 400; Real center_right_y = -40;

  if (x_coord >= center_left_x - 55 and x_coord <= center_left_x + 55 and y_coord >= center_left_y - 10 and y_coord <= center_left_y + 10)
  {
    //T1_o = 82.0524e6; //1% overstress
    T1_o = 81.6e6;
  }
  else if (x_coord >= center_right_x - 55 and x_coord <= center_right_x + 55 and y_coord >= center_right_y - 10 and y_coord <= center_right_y + 10)
  {
    T1_o = 81.6e6;
  }
  else{
    T1_o = 70e6;
  }
  
  return T1_o;

}