#include "SpatialDamageBreakageParameters.h"

registerMooseObject("farmsApp", SpatialDamageBreakageParameters);

InputParameters 
SpatialDamageBreakageParameters::validParams()
{
  InputParameters params = Function::validParams();
  params.addClassDescription("A function that defines spatial damage breakage parameters.");
  params.addRequiredParam<Real>("xmin", "The minimum x-coordinate.");
  params.addRequiredParam<Real>("xmax", "The maximum x-coordinate.");
  params.addRequiredParam<Real>("ymin", "The minimum y-coordinate.");
  params.addRequiredParam<Real>("ymax", "The maximum y-coordinate.");
  params.addRequiredParam<Real>("max_val", "The maximum value.");
  params.addRequiredParam<Real>("min_val", "The minimum value."); // Add this line
  params.addRequiredParam<Real>("scale", "The scale factor.");
  return params;
}

SpatialDamageBreakageParameters::SpatialDamageBreakageParameters(const InputParameters & parameters)
  : Function(parameters),
    _xmin(getParam<Real>("xmin")),
    _xmax(getParam<Real>("xmax")),
    _ymin(getParam<Real>("ymin")),
    _ymax(getParam<Real>("ymax")),
    _max_val(getParam<Real>("max_val")),
    _min_val(getParam<Real>("min_val")), // Add this line
    _scale(getParam<Real>("scale"))
{
}

Real
SpatialDamageBreakageParameters::value(Real /*t*/, const Point & p) const
{

 // Coordinates
  const Real x = p(0);
  const Real y = p(1);

  // 1) Compute how far (in Euclidian distance) we are from the box if we are outside it;
  //    inside the box => distance = 0
  Real dx = 0.0;
  if      (x < _xmin) dx = _xmin - x;
  else if (x > _xmax) dx = x - _xmax;

  Real dy = 0.0;
  if      (y < _ymin) dy = _ymin - y;
  else if (y > _ymax) dy = y - _ymax;

  // Euclidian distance outside the box
  Real dist = std::sqrt(dx * dx + dy * dy);

  // 2) Inside the box => dist=0 => exp(0) = 1.0
  //    Outside => decays exponentially with distance
  //    If you prefer a steeper or gentler transition, change the exponent!
  Real param_value = std::max(_min_val, _max_val * std::exp(-dist / _scale));

  return param_value;  

}