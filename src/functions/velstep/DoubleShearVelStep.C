#include "DoubleShearVelStep.h"
#include <string>

registerMooseObject("farmsApp", DoubleShearVelStep);

InputParameters
DoubleShearVelStep::validParams()
{
  InputParameters params = Function::validParams();

  // User-defined durations for each piece of the velocity function
  params.addParam<Real>("T0", 1.0, "Time duration for initial steady velocity segment");
  params.addParam<Real>("T1", 1.0, "Time duration for ramp from v_steady to v_high");
  params.addParam<Real>("T2", 1.0, "Time duration for ramp from v_high to v_low");

  // User-defined velocities
  params.addParam<Real>("v1", 1.0,  "Steady velocity in Segment 1");
  params.addParam<Real>("v2",   100.0, "Velocity after first ramp");
  params.addParam<Real>("v3",    10.0,  "Velocity after second ramp");

  return params;
}

DoubleShearVelStep::DoubleShearVelStep(const InputParameters & parameters)
  : Function(parameters),
    _T0(getParam<Real>("T0")),
    _T1(getParam<Real>("T1")),
    _T2(getParam<Real>("T2")),
    _v1(getParam<Real>("v1")),
    _v2(getParam<Real>("v2")),
    _v3(getParam<Real>("v3"))
{
}

Real
DoubleShearVelStep::value(Real t, const Point & /*p*/) const
{
  // We define piecewise segments in time:
  //
  // Segment 1:  0 <= t <= T0
  //    velocity = v_steady
  //    displacement d(t) = v_steady * t
  //
  // Segment 2:  T0 < t <= T0 + T1
  //    velocity ramps linearly from v_steady to v_high
  //    v(t) = v_steady + (v_high - v_steady)*((t - T0)/T1)
  //    d(t) = d(T0) + ∫(T0..t) v(τ) dτ
  //
  // Segment 3:  T0 + T1 < t <= T0 + T1 + T2
  //    velocity ramps linearly from v_high to v_low
  //    v(t) = v_high + (v_low - v_high)*((t - (T0 + T1))/T2)
  //    d(t) = d(T0 + T1) + ∫(T0+T1..t) v(τ) dτ
  //
  // Segment 4:  t >= T0 + T1 + T2
  //    velocity = v_low
  //    d(t) = d(T0 + T1 + T2) + v_low * (t - (T0 + T1 + T2))

  // For clarity, define time boundaries:
  const Real t1_end = _T0;         // End of Segment 1
  const Real t2_end = _T0 + _T1;   // End of Segment 2
  // const Real t3_end = _T0 + _T1 + _T2; // End of Segment 3

  Real disp = 0.0;

  // --- Segment 1: 0 <= t <= T0 ---
  if (t <= t1_end)
  {
    // d(t) = v_steady * t
    disp = _v1 * t;
  }
  // --- Segment 2: T0 < t <= T0 + T1 ---
  else if (t <= t2_end)
  {
    Real disp_T0 = _v1 * t1_end;

    Real x = t - t1_end;
    Real disp_segment2 = _v2 * x;

    disp = disp_T0 + disp_segment2;
  }
  // --- Segment 3: T0 + T1 < t <= T0 + T1 + T2 ---
  else
  {
    Real disp_T0 = _v1 * t1_end;
    Real disp_T1 = _v2 * ( t2_end - t1_end );

    Real disp_T1_end = disp_T0 + disp_T1;

    Real x = t - t2_end;
    Real disp_segment3 = _v3 * x;

    disp = disp_T1_end + disp_segment3;
  }

  return disp;
}