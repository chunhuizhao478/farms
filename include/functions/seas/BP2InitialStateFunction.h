//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Function.h"

/**
 * Computes initial state variable θ(z,0) for BP2 benchmark using Eq. 12:
 * θ(z,0) = (Dc/V0) * exp{(a(z)/b) * ln[2V0/Vinit * sinh((τ⁰ - η*Vinit)/(a(z)*σn))] - f0/b}
 *
 * where a(z) is depth-dependent:
 *   a = a0 for z < H
 *   a = a0 + (amax-a0)*(z-H)/h for H <= z < H+h
 *   a = amax for z >= H+h
 */
class BP2InitialStateFunction : public Function
{
public:
  static InputParameters validParams();

  BP2InitialStateFunction(const InputParameters & parameters);

  virtual Real value(Real t, const Point & p) const override;

protected:
  /// Compute a(z) based on depth
  Real computeA(Real z) const;

  /// Critical slip distance
  const Real _Dc;
  /// Reference slip rate
  const Real _V0;
  /// Initial slip rate
  const Real _Vinit;
  /// Rate-state parameter b
  const Real _b;
  /// Reference friction coefficient
  const Real _f0;
  /// Pre-stress
  const Real _tau0;
  /// Normal stress
  const Real _sigma_n;
  /// Radiation damping coefficient
  const Real _eta;
  /// VW region a value
  const Real _a0;
  /// VS region a value
  const Real _amax;
  /// Depth of VW region
  const Real _H;
  /// Width of transition zone
  const Real _h;
  /// Coordinate direction for depth (default: 1 = y)
  const unsigned int _depth_dir;
};
