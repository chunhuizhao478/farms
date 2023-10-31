/*
Define Function for Initial Shear Stress Hetergeneoity in XY Direction 
Will be Passed into AuxVariables
Chunhui Zhao

cluster run stepover_3
*/

#include "InitialStressXYstepover3.h"

registerMooseObject("farmsApp", InitialStressXYstepover3);

InputParameters
InitialStressXYstepover3::validParams()
{
  InputParameters params = Function::validParams();
  return params;
}

InitialStressXYstepover3::InitialStressXYstepover3(const InputParameters & parameters)
  : Function(parameters)
{
}

Real
InitialStressXYstepover3::value(Real /*t*/, const Point & p) const
{
  Real x_coord = p(0);
  Real y_coord = p(1);
  Real regiondx = 200; //mu_d=0.49 //Real regiondx = 100; //mu_d=0.1
  Real T1_o = 0.0;
  if ((x_coord<=(-4000+regiondx))&&(x_coord>=(-4000-regiondx))&&(y_coord<=(0.0+regiondx))&&(y_coord>=(0.0-regiondx)))
  {
    //T1_o = 81.6e6;
    //T1_o = 84.5e6; //0.04
    T1_o = 82.0524e6; //0.01
  }
  else
  {
    T1_o = 70e6;
  }
  return T1_o;


}