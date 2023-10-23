/*
Define Function for Initial Shear Stress Hetergeneoity in XY Direction 
Will be Passed into AuxVariables
Chunhui Zhao

cluster run stepover_3
*/

#include "InitialStressXYfinitelen_multiactivation.h"

registerMooseObject("farmsApp", InitialStressXYfinitelen_multiactivation);

InputParameters
InitialStressXYfinitelen_multiactivation::validParams()
{
  InputParameters params = Function::validParams();
  return params;
}

InitialStressXYfinitelen_multiactivation::InitialStressXYfinitelen_multiactivation(const InputParameters & parameters)
  : Function(parameters)
{
}

Real
InitialStressXYfinitelen_multiactivation::value(Real t, const Point & p) const
{
  Real x_coord = p(0);
  Real y_coord = p(1);
  Real T1_o = 0.0;

  if ( t >= 0.0 && t <= 2.0 ){
    if ((x_coord<=(-4000+100))&&(x_coord>=(-4000-100))&&(y_coord<=(0.0+100))&&(y_coord>=(0.0-100)))
    {
      T1_o = 81.6e6;
    }
    else
    {
      T1_o = 70e6;
    }
  }
  else if ( t > 2.0 && t <= 4.0 ){
    if ((x_coord<=(-4000+100))&&(x_coord>=(-4000-100))&&(y_coord<=(0.0+100))&&(y_coord>=(0.0-100)))
    {
      T1_o = 81.6e6;
    }
    else if ((x_coord<=(-2500+100))&&(x_coord>=(-2500-100))&&(y_coord<=(0.0+100))&&(y_coord>=(0.0-100)))
    {
      T1_o = 81.6e6;
    }
    else
    {
      T1_o = 70e6;
    }
  }
  else{
    if ((x_coord<=(-4000+100))&&(x_coord>=(-4000-100))&&(y_coord<=(0.0+100))&&(y_coord>=(0.0-100)))
    {
      T1_o = 81.6e6;
    }
    else if ((x_coord<=(-3000+100))&&(x_coord>=(-3000-100))&&(y_coord<=(0.0+100))&&(y_coord>=(0.0-100)))
    {
      T1_o = 81.6e6;
    }
    else if ((x_coord<=(-2000+100))&&(x_coord>=(-2000-100))&&(y_coord<=(0.0+100))&&(y_coord>=(0.0-100)))
    {
      T1_o = 81.6e6;
    }
    else
    {
      T1_o = 70e6;
    }
  }
  return T1_o;

  // //use Gaussian
  // Real T1_o = 0.0;

  // Real x_coord = p(0); //along the strike direction
  // Real y_coord = p(1); //along the normal direction

  // //Parameter
  // Real T1_perturb = 11.6e6; //pertubation (Pa) 81.6e6-70e6
  // Real R = 500; //(m)

  // Real r = sqrt((x_coord+2000)*(x_coord+2000)); //2D

  // Real Val_F = 0.0;
  // if (r < R){
  //   Val_F = exp((r*r)/(r*r-R*R));
  // }
  // else{
  //   Val_F = 0.0;
  // }

  // Real T1_perturb_o = T1_perturb * Val_F;

  // if ((y_coord<=(0.0+100))&&(y_coord>=(0.0-100))){
  //   T1_o = 70e6 + T1_perturb_o;
  // }
  // else{
  //   T1_o = 70e6; 
  // }

  // return T1_o;

}