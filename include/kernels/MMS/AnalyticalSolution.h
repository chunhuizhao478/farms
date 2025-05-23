#pragma once
#include "Function.h"

/**
 * Provides analytical solution for MMS verification
 */
class AnalyticalSolution : public Function
{
public:
  static InputParameters validParams();
  AnalyticalSolution(const InputParameters & parameters);

  virtual Real value(Real t, const Point & p) const override;

protected:
  // Parameters
  const Real _amplitude;
  const Real _sigma;
  const Real _x0;
  const Real _y0;
  const Real _t0;
  const Real _t_width;
  const Real _omega;
  const unsigned int _component;
};