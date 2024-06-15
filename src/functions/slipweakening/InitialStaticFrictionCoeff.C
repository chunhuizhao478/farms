/*
Define Function for Initial Static Friction Coefficient for TPV24 benchmark
*/

#include "InitialStaticFrictionCoeff.h"

registerMooseObject("farmsApp", InitialStaticFrictionCoeff);

InputParameters
InitialStaticFrictionCoeff::validParams()
{
  InputParameters params = Function::validParams();
  return params;
}

InitialStaticFrictionCoeff::InitialStaticFrictionCoeff(const InputParameters & parameters)
  : Function(parameters)
{
}

Real
InitialStaticFrictionCoeff::value(Real /*t*/, const Point & p) const
{

  Real x_coord = p(0); //along the x direction
  // Real y_coord = p(1); //along the y direction
  Real z_coord = p(2); //along the z direction

  //TPV205-3D
  Real mu_s = 0.0;
  if ((x_coord<=(0.0+15e3))&&(x_coord>=(0.0-15e3)) && (z_coord<=(0.0+0.0))&&(z_coord>=(0.0-15e3)))
  {
	  mu_s = 0.677;
  }  
  else
  {
    mu_s = 10000.0;
  }

  return mu_s;

}