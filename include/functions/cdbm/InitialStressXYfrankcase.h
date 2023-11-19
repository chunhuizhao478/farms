/*
Define Function for Initial Shear Stress Hetergeneoity in XY Direction 
frank case
Chunhui Zhao
*/

#pragma once

#include "Function.h"

class InitialStressXYfrankcase : public Function
{
public:
  InitialStressXYfrankcase(const InputParameters & parameters);

  static InputParameters validParams();

  using Function::value;
  virtual Real value(Real t, const Point & p) const override;

};