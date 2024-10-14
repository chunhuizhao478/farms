#include "ForcedRuptureTimeTPV243D.h"

registerMooseObject("farmsApp", ForcedRuptureTimeTPV243D);

InputParameters
ForcedRuptureTimeTPV243D::validParams()
{
  InputParameters params = Function::validParams();
  params.addRequiredParam<Real>("loc_x", "nucleation point x coordinate");
  params.addRequiredParam<Real>("loc_y", "nucleation point y coordinate");
  params.addRequiredParam<Real>("loc_z", "nucleation point z coordinate"); 
  params.addRequiredParam<Real>("r_crit", "critical distance to hypocenter");
  params.addRequiredParam<Real>("Vs", "shear wave speed");   
  return params;
}

ForcedRuptureTimeTPV243D::ForcedRuptureTimeTPV243D(const InputParameters & parameters)
  : Function(parameters),
  _loc_x(getParam<Real>("loc_x")),
  _loc_y(getParam<Real>("loc_y")),
  _loc_z(getParam<Real>("loc_z")),
  _r_crit(getParam<Real>("r_crit")),
  _Vs(getParam<Real>("Vs"))   
{
}

Real
ForcedRuptureTimeTPV243D::value(Real /*t*/, const Point & p) const
{
  
  //the coordinate follows benchmark
  Real x_coord = p(0); //along the strike direction
  Real y_coord = p(1); //along the dip direction
  Real z_coord = p(2); //along the normal direction

  Real r = sqrt( (x_coord - _loc_x) * (x_coord - _loc_x) + (y_coord - _loc_y) * (y_coord - _loc_y) + (z_coord - _loc_z) * (z_coord - _loc_z) );
  
  Real T = 0.0;
  if ( r < _r_crit ){
    T = r / ( 0.7 * _Vs ) + ( 0.081 * _r_crit ) / ( 0.7 * _Vs ) * ( 1 / (1 - (r/_r_crit)*(r/_r_crit)) - 1 );
  }
  else{
    T = 1e9;
  }

  return T;

}