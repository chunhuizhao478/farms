#pragma once

#include "Function.h"

class InitialShearStress3D : public Function
{
public:
  InitialShearStress3D(const InputParameters & parameters);

  static InputParameters validParams();

  using Function::value;
  virtual Real value(Real t, const Point & p) const override;

};