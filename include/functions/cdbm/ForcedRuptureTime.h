/*
Define Forced Rupture Time T (TPV24)
Gaussian Distribution
Chunhui Zhao
*/

#pragma once

#include "Function.h"

class ForcedRuptureTime : public Function
{
public:
  ForcedRuptureTime(const InputParameters & parameters);

  static InputParameters validParams();

  using Function::value;
  virtual Real value(Real t, const Point & p) const override;

  Real _loc_x;
  Real _loc_y;
  Real _loc_z;
  Real _r_crit;
  Real _Vs;

};