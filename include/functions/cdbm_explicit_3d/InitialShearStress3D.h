#pragma once

#include "Function.h"

class InitialShearStress3D : public Function
{
public:
  InitialShearStress3D(const InputParameters & parameters);

  static InputParameters validParams();

  using Function::value;
  virtual Real value(Real t, const Point & p) const override;

  /// nucleation center (x,y,z)
  std::vector<Real> _nucl_center;

  /// sigma (exponential decay)
  Real _sigma;  

  /// sigma (exponential decay)
  Real _min_val;

  Real _max_val;  

};