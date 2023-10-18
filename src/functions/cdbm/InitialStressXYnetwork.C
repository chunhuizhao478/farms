/*
Define Function for Initial Shear Stress Hetergeneoity in XY Direction 
Will be Passed into AuxVariables
Chunhui Zhao
*/

#include "InitialStressXYnetwork.h"

registerMooseObject("farmsApp", InitialStressXYnetwork);

InputParameters
InitialStressXYnetwork::validParams()
{
  InputParameters params = Function::validParams();
  return params;
}

InitialStressXYnetwork::InitialStressXYnetwork(const InputParameters & parameters)
  : Function(parameters)
{
}

Real
InitialStressXYnetwork::value(Real /*t*/, const Point & p) const
{
  // Real x_coord = p(0);
  // Real y_coord = p(1);
  // double T1_o = 0.0;
  // if ((x_coord<=(1000.0+2e3))&&(x_coord>=(1000.0-2e3))&&(y_coord<=(-1000.0+2e3))&&(y_coord>=(-1000.0-2e3)))
  // {
  //     T1_o = 81.6e6;
  // }
  // else
  // {
  //     T1_o = 70.0e6;
  // }
  // return T1_o;

  //not center
  // Real x_coord = p(0);
  // Real y_coord = p(1);
  // double T1_o = 0.0;
  // if ((x_coord<=(-300+1.5e3))&&(x_coord>=(-300+1.5e3))&&(y_coord<=(1000.0+400))&&(y_coord>=(1000.0-400)))
  // {
  //     T1_o = 81.6e6;
  // }
  // else
  // {
  //     T1_o = 70.0e6;
  // }
  // return T1_o;

  // //use Gaussian 
  // Real T1_o = 0.0;

  Real x_coord = p(0); //along the strike direction
  Real y_coord = p(1); //along the normal direction

  // //Parameter
  // Real T1_perturb = 11.6e6; //pertubation (Pa) 81.6e6-70e6
  // Real R = 2000; //(m)
  // Real center_x = 0;
  // Real center_y = 0;

  // Real r = sqrt( (x_coord - center_x) * (x_coord - center_x) + (y_coord - center_y) * (y_coord - center_y) ); //2D

  // Real Val_F = 0.0;
  // if (r < R){
  //   Val_F = exp((r*r)/(r*r-R*R));
  // }
  // else{
  //   Val_F = 0.0;
  // }

  // Real T1_perturb_o = T1_perturb * Val_F;

  // T1_o = 70e6 + T1_perturb_o;

  // return T1_o;

  Real T1_o = 0;

  if (x_coord >= -50 -100 and x_coord <= -50 + 100 and y_coord >= 110 -100 and y_coord <= 110 + 100)
  {
    T1_o = 81.6e6;
  }
  else{
    T1_o = 70e6;
  }
  
  return T1_o;

}