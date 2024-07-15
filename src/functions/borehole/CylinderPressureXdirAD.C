#include "CylinderPressureXdirAD.h"

registerMooseObject("farmsApp", CylinderPressureXdirAD);

InputParameters
CylinderPressureXdirAD::validParams()
{
  InputParameters params = Function::validParams();
  params.addRequiredParam<Real>("value","pressure value");
  return params;
}

CylinderPressureXdirAD::CylinderPressureXdirAD(const InputParameters & parameters)
  : Function(parameters),
  _value(getParam<Real>("value"))
{
}

Real
CylinderPressureXdirAD::value(Real /*t*/, const Point & p) const
{


  Real x_coord = p(0); 
  Real y_coord = p(1); 

  Real pressure_x = 0.0;

  Real mag = sqrt(pow(x_coord,2)+pow(y_coord,2));

  pressure_x = _value * abs( x_coord / mag );
  
  return pressure_x;

}