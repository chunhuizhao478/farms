#include "SpatialDamageBreakageParameters.h"

registerMooseObject("farmsApp", SpatialDamageBreakageParameters);

InputParameters 
SpatialDamageBreakageParameters::validParams()
{
  InputParameters params = Function::validParams();
  params.addClassDescription("A function that defines spatial damage breakage parameters.");
  params.addRequiredParam<Real>("W", "The width where the maximum value is defined.");
  params.addRequiredParam<Real>("w", "The transition width the maximum value reduces to the minimum value.");
  params.addRequiredParam<Real>("max_val", "The maximum value.");
  params.addRequiredParam<Real>("min_val", "The minimum value."); // Add this line
  return params;
}

SpatialDamageBreakageParameters::SpatialDamageBreakageParameters(const InputParameters & parameters)
  : Function(parameters),
    _W(getParam<Real>("W")),
    _w(getParam<Real>("w")),
    _max_val(getParam<Real>("max_val")),
    _min_val(getParam<Real>("min_val")) // Add this line
{
}

Real
SpatialDamageBreakageParameters::value(Real /*t*/, const Point & p) const
{

  Real val = 0.0;

  // Coordinates
  const Real x = p(0);
  const Real y = p(1);

  if ( std::abs(x) < _W ){
    val = _max_val;
  }
  else if ( std::abs(x) > _W + _w ){
    val = _min_val;
  }
  else if ( _w > 0 ){
    val = _max_val - (_max_val - _min_val) * (std::abs(x) - _W) / _w;
  }

  return val;

}