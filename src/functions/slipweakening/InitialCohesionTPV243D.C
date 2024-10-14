#include "InitialCohesionTPV243D.h"

registerMooseObject("farmsApp", InitialCohesionTPV243D);

InputParameters
InitialCohesionTPV243D::validParams()
{
  InputParameters params = Function::validParams();
  return params;
}

InitialCohesionTPV243D::InitialCohesionTPV243D(const InputParameters & parameters)
  : Function(parameters)
{
}

Real
InitialCohesionTPV243D::value(Real /*t*/, const Point & p) const
{
  
  //the coordinate follows benchmark
  Real y_coord = p(1); //along the dip direction
  Real Co = 0;

  if ( abs(y_coord) <= 4000 ){
    Co = 0.30 * 1e6 + ( 0.000675 * 1e6 ) * ( 4000 - abs(y_coord) );
  }
  else{
    Co = 0.30 * 1e6;
  }

  return Co;

}