/*
Define Function for Initial Static Friction Coefficient for benchmark
*/

#pragma once

#include "Function.h"

class InitialDipStressCDBM : public Function
{
public:
  InitialDipStressCDBM(const InputParameters & parameters);

  static InputParameters validParams();

  using Function::value;
  virtual Real value(Real t, const Point & p) const override;

  std::string _benchmark;

};