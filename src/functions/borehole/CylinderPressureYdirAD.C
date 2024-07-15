#include "CylinderPressureYdirAD.h"

registerMooseObject("farmsApp", CylinderPressureYdirAD);

InputParameters
CylinderPressureYdirAD::validParams()
{
  InputParameters params = Function::validParams();
  params.addRequiredParam<Real>("value","pressure value");
  return params;
}

CylinderPressureYdirAD::CylinderPressureYdirAD(const InputParameters & parameters)
  : Function(parameters),
  _value(getParam<Real>("value"))
{
}

Real
CylinderPressureYdirAD::value(Real /*t*/, const Point & p) const
{


  Real x_coord = p(0); 
  Real y_coord = p(1); 

  Real pressure_y = 0.0;

  Real mag = sqrt(pow(x_coord,2)+pow(y_coord,2));

  pressure_y = _value * abs( y_coord / mag );
  
  return pressure_y;

}