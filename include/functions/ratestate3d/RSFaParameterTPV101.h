/*
Define Function for Spatially Variable RSF 'a' Parameter
Problem-Specific: TPV101-3D

According to SCEC TPV101 benchmark:
a(x,y) = a_0 + Δa(x,y)
Δa(x,y) = Δa_0 * [1 - B(x;W,w) * B(y-y_0; W/2, w)]

where B is the smooth boxcar function.

Coordinate system in this implementation:
- x: along-strike direction
- y: fault-normal direction (fault at y=0)
- z: depth direction (z=0 at surface, z<0 is depth)

Mapping to benchmark coordinates:
- Implementation x → Benchmark x (along-strike)
- Implementation |z| → Benchmark y (depth)
*/

#pragma once

#include "Function.h"

class RSFaParameterTPV101 : public Function
{
public:
  RSFaParameterTPV101(const InputParameters & parameters);

  static InputParameters validParams();

  using Function::value;
  virtual Real value(Real t, const Point & p) const override;

protected:
  /// Smooth boxcar function B(x; W, w)
  Real boxcarB(Real x, Real W, Real w) const;

  /// Base value of a in velocity-weakening region
  const Real _a_0;

  /// Maximum increase in a for velocity-strengthening region
  const Real _delta_a_0;

  /// Half-width of velocity-weakening region (m)
  const Real _W;

  /// Transition layer width (m)
  const Real _w;

  /// Depth of hypocenter / center of VW region in depth direction (m)
  const Real _y_0;
};
