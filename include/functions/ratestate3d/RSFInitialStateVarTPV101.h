/*
Define Function for Spatially Variable Initial State Variable
Problem-Specific: TPV101-3D

According to SCEC TPV101 benchmark (Eq. 6):
θ_ini(x,y) = (L/V_0) * exp[(a*ln(2*sinh(τ_ini/(a*σ_ini))) - f_0 - a(x,y)*ln(V_ini/V_0)) / b]

The initial state variable must vary spatially to maintain uniform initial
shear stress (τ_ini = 75 MPa) given the spatially variable 'a' parameter.

Coordinate system in this implementation:
- x: along-strike direction
- y: fault-normal direction (fault at y=0)
- z: depth direction (z=0 at surface, z<0 is depth)
*/

#pragma once

#include "Function.h"

class RSFInitialStateVarTPV101 : public Function
{
public:
  RSFInitialStateVarTPV101(const InputParameters & parameters);

  static InputParameters validParams();

  using Function::value;
  virtual Real value(Real t, const Point & p) const override;

protected:
  /// Smooth boxcar function B(x; W, w)
  Real boxcarB(Real x, Real W, Real w) const;

  /// Compute local a(x,y) value
  Real computeLocalA(const Point & p) const;

  /// RSF parameters
  const Real _f_0;      // Reference friction coefficient
  const Real _V_0;      // Reference slip velocity (m/s)
  const Real _a_0;      // Base value of a in VW region
  const Real _b;        // Evolution effect parameter
  const Real _L;        // Characteristic slip distance (m)
  const Real _delta_a_0; // Maximum increase in a for VS region

  /// Initial conditions
  const Real _tau_ini;  // Initial shear stress (Pa)
  const Real _sigma_ini; // Initial normal stress (Pa)
  const Real _V_ini;    // Initial slip velocity (m/s)

  /// Geometry parameters
  const Real _W;        // Half-width of VW region (m)
  const Real _w;        // Transition layer width (m)
  const Real _y_0;      // Depth of center of VW region (m)
};
