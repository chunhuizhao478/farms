/*
Define Function for Initial Static Friction Coefficient for TPV24 benchmark
*/

#pragma once

#include "Function.h"

class InitialStaticFrictionCoeffTPV243D : public Function
{
public:
  InitialStaticFrictionCoeffTPV243D(const InputParameters & parameters);

  static InputParameters validParams();

  using Function::value;
  virtual Real value(Real t, const Point & p) const override;

};