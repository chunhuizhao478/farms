/*
Define Function for Initial Shear Stress Hetergeneoity in XY Direction 
Gaussian Distribution
Chunhui Zhao
*/

#pragma once

#include "Function.h"

class InitialAlphaAD : public Function
{
public:
  InitialAlphaAD(const InputParameters & parameters);

  static InputParameters validParams();

  using Function::value;
  virtual Real value(Real t, const Point & p) const override;

};