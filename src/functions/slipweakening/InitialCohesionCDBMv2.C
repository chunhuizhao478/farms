#include "InitialCohesionCDBMv2.h"

registerMooseObject("farmsApp", InitialCohesionCDBMv2);

InputParameters
InitialCohesionCDBMv2::validParams()
{
  InputParameters params = Function::validParams();
  params.addClassDescription("Initial cohesion function for CDBMv2 benchmark with depth-dependent "
                             "cohesion and end-zone strengthening.");
  params.addRequiredParam<Real>("depth", "The depth at which the cohesion starts to decrease, in meters.");
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
  // Fault geometry (hardcoded)
  const Real fault_half_length = 22500.0;  // Fault extends from -22500 to +22500 m
  const Real taper_zone_width = 2500.0;    // Last 2500 m on each end
  const Real target_cohesion_increase = 1.0e8; // 100 MPa increase in taper zones (in Pa)
  
  // Depth-dependent cohesion (original implementation)
  Real z_coord = p(2); // Along the dip direction (depth)
  Real Co = 0.0;
  
  if (std::abs(z_coord) <= _depth)
  {
    Co = _min_cohesion * 1e6 + (_slope * 1e6) * (_depth - std::abs(z_coord));
  }
  else
  {
    Co = _min_cohesion * 1e6;
  }
  
  // Lateral coordinate along fault strike (x-direction)
  Real x_coord = p(0);
  
  // Calculate end-zone strengthening
  Real cohesion_modifier = 0.0;
  
  // Right end taper zone (20000 < x <= 22500)
  if (x_coord > (fault_half_length - taper_zone_width) && x_coord <= fault_half_length)
  {
    Real dist_into_taper = x_coord - (fault_half_length - taper_zone_width);
    Real normalized_dist = dist_into_taper / taper_zone_width; // 0 to 1
    cohesion_modifier = target_cohesion_increase * normalized_dist;
  }
  // Left end taper zone (-22500 <= x < -20000)
  else if (x_coord < -(fault_half_length - taper_zone_width) && x_coord >= -fault_half_length)
  {
    Real dist_into_taper = std::abs(x_coord) - (fault_half_length - taper_zone_width);
    Real normalized_dist = dist_into_taper / taper_zone_width; // 0 to 1
    cohesion_modifier = target_cohesion_increase * normalized_dist;
  }
  // Middle section of fault: no modification
  
  return Co + cohesion_modifier;
}