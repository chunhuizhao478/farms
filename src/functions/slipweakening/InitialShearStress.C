/*
Define Function for Initial Shear Stress for benchmark
*/

#include "InitialShearStress.h"

registerMooseObject("farmsApp", InitialShearStress);

InputParameters
InitialShearStress::validParams()
{
  InputParameters params = Function::validParams();
  return params;
}

InitialShearStress::InitialShearStress(const InputParameters & parameters)
  : Function(parameters)
{
}

Real
InitialShearStress::value(Real /*t*/, const Point & p) const
{

  Real x_coord = p(0); //along the x direction
  // Real y_coord = p(1); //along the y direction
  Real z_coord = p(2); //along the z direction

  Real T1_o = 0.0;
  if ((x_coord<=(0.0+1.5e3))&&(x_coord>=(0.0-1.5e3))&& (z_coord<=(-7.5e3+1.5e3))&&(z_coord>=(-7.5e3-1.5e3)))
  {
      T1_o = 81.6e6;
  }
  else if ((x_coord<=(-7.5e3+1.5e3))&&(x_coord>=(-7.5e3-1.5e3)) && (z_coord<=(-7.5e3+1.5e3))&&(z_coord>=(-7.5e3-1.5e3)))
  {
      T1_o = 78.0e6;
  }
  else if ((x_coord<=(7.5e3+1.5e3))&&(x_coord>=(7.5e3-1.5e3)) && (z_coord<=(-7.5e3+1.5e3))&&(z_coord>=(-7.5e3-1.5e3)))
  {
      T1_o = 62.0e6;
  }
  else
  {
      T1_o = 70.0e6;
  }
  return T1_o;

}