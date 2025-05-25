#pragma once

#include "Function.h"

class InitialDamageCycleSim3DPlaneFunction : public Function
{
public:
  InitialDamageCycleSim3DPlaneFunction(const InputParameters & parameters);

  static InputParameters validParams();

  using Function::value;
  virtual Real value(Real t, const Point & p) const override;

  /// Length of the fault in the x-direction
  Real _len_of_fault_strike;

  /// Length of the fault in the z-direction (the 'dip' extent)
  Real _len_of_fault_dip;

  /// sigma (exponential decay)
  Real _sigma;

  /// Peak value of the initial damage at the plane
  Real _peak_val;

  /// nucleation center (x,y,z)
  std::vector<Real> _nucl_center;

};