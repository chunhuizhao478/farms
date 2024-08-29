/*
Define Function for Initial Static Friction Coefficient for TPV24 benchmark
*/

#include "InitialStaticFrictionCoeffTPV243D.h"

registerMooseObject("farmsApp", InitialStaticFrictionCoeffTPV243D);

InputParameters
InitialStaticFrictionCoeffTPV243D::validParams()
{
  InputParameters params = Function::validParams();
  return params;
}

InitialStaticFrictionCoeffTPV243D::InitialStaticFrictionCoeffTPV243D(const InputParameters & parameters)
  : Function(parameters)
{
}

Real
InitialStaticFrictionCoeffTPV243D::value(Real /*t*/, const Point & /*p*/) const
{

  // Real x_coord = p(0); //along the x direction
  // Real y_coord = p(1); //along the y direction
  // Real z_coord = p(2); //along the z direction

  Real mu_s = 0;

  //Parameter
  //if ( z_coord < -15000 || x_coord < -16000 || x_coord > 12000 || ( x_coord > 10392.5 && y_coord < -6000 ) ){
  //  mu_s = 1e12;
  //}
  //else{
    mu_s = 0.18;
  //}

  return mu_s;

}