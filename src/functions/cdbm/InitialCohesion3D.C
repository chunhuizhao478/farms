#include "InitialCohesion3D.h"

registerMooseObject("farmsApp", InitialCohesion3D);

InputParameters
InitialCohesion3D::validParams()
{
  InputParameters params = Function::validParams();
  return params;
}

InitialCohesion3D::InitialCohesion3D(const InputParameters & parameters)
  : Function(parameters)
{
}

Real
InitialCohesion3D::value(Real /*t*/, const Point & p) const
{
  
  //the coordinate follows benchmark
  Real z_coord = p(2); //along the slip direction
  Real Co = 0;

  if ( abs(z_coord) <= 4000 ){
    Co = 0.30 * 1e6 + ( 0.000675 * 1e6 ) * ( 4000 - abs(z_coord) );
  }
  else{
    Co = 0.30 * 1e6;
  }

  return Co;

}