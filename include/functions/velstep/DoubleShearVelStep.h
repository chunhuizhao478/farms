#pragma once

#include "Function.h"

/**
 * A Function that provides a piecewise velocity profile:
 *
 * Segment 1: 0 <= t <= T0
 *     velocity = v_steady
 *
 * Segment 2: T0 < t <= T0 + T1
 *     velocity ramps linearly from v_steady to v_high
 *
 * Segment 3: T0 + T1 < t <= T0 + T1 + T2
 *     velocity ramps linearly from v_high to v_low
 *
 * Segment 4: t >= T0 + T1 + T2
 *     velocity = v_low
 */
class DoubleShearVelStep : public Function
{
public:
  static InputParameters validParams();
  DoubleShearVelStep(const InputParameters & parameters);

  virtual Real value(Real t, const Point & p) const override;

protected:

  // Durations for each piecewise segment
  Real _T0;  // Duration of steady velocity
  Real _T1;  // Duration of ramp from v_steady to v_high
  Real _T2;  // Duration of ramp from v_high to v_low

  // Velocities
  Real _v1; // Steady velocity in Segment 1
  Real _v2;   // Velocity at the end of ramp-up
  Real _v3;    // Velocity after ramp-down
};