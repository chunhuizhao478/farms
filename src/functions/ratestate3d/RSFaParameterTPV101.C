/*
Define Function for Spatially Variable RSF 'a' Parameter
Problem-Specific: TPV101-3D

According to SCEC TPV101 benchmark (Eq. 4-5):
a(x,y) = a_0 + Δa(x,y)
Δa(x,y) = Δa_0 * [1 - B(x;W,w) * B(y-y_0; W/2, w)]

B(x;W,w) = {
  1,                                           |x| <= W
  0.5*[1 + tanh(w/(|x|-W-w) + w/(|x|-W))],    W < |x| < W+w
  0,                                           |x| >= W+w
}
*/

#include "RSFaParameterTPV101.h"
#include <cmath>

registerMooseObject("farmsApp", RSFaParameterTPV101);

InputParameters
RSFaParameterTPV101::validParams()
{
  InputParameters params = Function::validParams();
  params.addClassDescription("Spatially variable RSF 'a' parameter for TPV101 benchmark.");
  params.addParam<Real>("a_0", 0.008, "Base value of a in velocity-weakening region");
  params.addParam<Real>("delta_a_0", 0.008, "Maximum increase in a for velocity-strengthening");
  params.addParam<Real>("W", 15000.0, "Half-width of velocity-weakening region (m)");
  params.addParam<Real>("w", 3000.0, "Transition layer width (m)");
  params.addParam<Real>("y_0", 7500.0, "Depth of center of VW region (m), positive value");
  return params;
}

RSFaParameterTPV101::RSFaParameterTPV101(const InputParameters & parameters)
  : Function(parameters),
    _a_0(getParam<Real>("a_0")),
    _delta_a_0(getParam<Real>("delta_a_0")),
    _W(getParam<Real>("W")),
    _w(getParam<Real>("w")),
    _y_0(getParam<Real>("y_0"))
{
}

Real
RSFaParameterTPV101::boxcarB(Real x, Real W, Real w) const
{
  Real abs_x = std::abs(x);

  if (abs_x <= W)
  {
    return 1.0;
  }
  else if (abs_x >= W + w)
  {
    return 0.0;
  }
  else
  {
    // Transition region: W < |x| < W + w
    // B = 0.5 * [1 + tanh(w/(|x|-W-w) + w/(|x|-W))]
    Real term1 = w / (abs_x - W - w);  // This is negative since |x| < W+w
    Real term2 = w / (abs_x - W);      // This is positive since |x| > W
    return 0.5 * (1.0 + std::tanh(term1 + term2));
  }
}

Real
RSFaParameterTPV101::value(Real /*t*/, const Point & p) const
{
  // Get coordinates
  // x: along-strike direction
  // z: depth direction (z <= 0, with z=0 at surface)
  Real x_coord = p(0);  // along-strike
  Real z_coord = p(2);  // depth (negative values)

  // Convert z to depth (positive value for benchmark convention)
  // In benchmark: y is depth with y=0 at surface, y>0 going down
  // In implementation: z is depth with z=0 at surface, z<0 going down
  Real depth = std::abs(z_coord);

  // Compute B functions
  // B_x = B(x; W, w) for along-strike direction
  Real B_x = boxcarB(x_coord, _W, _w);

  // B_y = B(y - y_0; W/2, w) for depth direction
  // y_0 = 7.5 km is the center of the VW region in depth
  // W/2 = 7.5 km is the half-width in depth direction
  Real B_y = boxcarB(depth - _y_0, _W / 2.0, _w);

  // Compute Δa(x,y) = Δa_0 * [1 - B_x * B_y]
  Real delta_a = _delta_a_0 * (1.0 - B_x * B_y);

  // Return a(x,y) = a_0 + Δa(x,y)
  return _a_0 + delta_a;
}
