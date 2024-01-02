/*
Define Function for Initial Shear Stress Hetergeneoity in XY Direction 
Will be Passed directly to "ComputeEigenstrainFromInitialStress" as variables
Chunhui Zhao
*/

#pragma once

#include "Function.h"

class InitialStressXYPressureBorehole_x : public Function
{
public:
  InitialStressXYPressureBorehole_x(const InputParameters & parameters);

  static InputParameters validParams();

  using Function::value;
  virtual Real value(Real t, const Point & p) const override;

};