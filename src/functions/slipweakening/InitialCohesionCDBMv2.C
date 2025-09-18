#include "InitialCohesionCDBMv2.h"

registerMooseObject("farmsApp", InitialCohesionCDBMv2);

InputParameters
InitialCohesionCDBMv2::validParams()
{
  InputParameters params = Function::validParams();
  params.addClassDescription("Initial cohesion function for CDBMv2 benchmark, defined by a piecewise linear function based on depth.");
  params.addRequiredParam<Real>("depth", "The depth at which the cohesion starts to decrease within, in meters.");
  params.addRequiredParam<Real>("slope", "The slope of the linear decrease in cohesion with depth, in MPa per meter.");
  params.addRequiredParam<Real>("min_cohesion", "The minimum cohesion value at the maximum depth, in MPa.");
  return params;
}

InitialCohesionCDBMv2::InitialCohesionCDBMv2(const InputParameters & parameters)
  : Function(parameters),
    _depth(getParam<Real>("depth")),
    _slope(getParam<Real>("slope")),
    _min_cohesion(getParam<Real>("min_cohesion"))
{
}

Real
InitialCohesionCDBMv2::value(Real /*t*/, const Point & p) const
{
  
  //the coordinate follows benchmark
  Real z_coord = p(2); //along the dip direction
  Real Co = 0;

  if ( abs(z_coord) <= _depth ){
    Co = _min_cohesion * 1e6 + ( _slope * 1e6 ) * ( _depth - abs(z_coord) );
  }
  else{
    Co = _min_cohesion * 1e6;
  }

  return Co;

}