#pragma once

#include "Function.h"

class SpatialDamageBreakageParameters : public Function
{
public:
  static InputParameters validParams();
  SpatialDamageBreakageParameters(const InputParameters & parameters);

  virtual Real value(Real t, const Point & p) const override;

private:
  const Real _xmin;
  const Real _xmax;
  const Real _ymin;
  const Real _ymax;
  const Real _max_val;
  const Real _min_val; // Add this line
  const Real _scale;
};