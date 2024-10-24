/*
Define Function for Initial Shear Stress Hetergeneoity in XY Direction 
Will be Passed directly to "ComputeEigenstrainFromInitialStress" as variables
Chunhui Zhao
*/

#pragma once

#include "Function.h"

class InitialStressXYnetwork : public Function
{
public:
  InitialStressXYnetwork(const InputParameters & parameters);

  static InputParameters validParams();

  using Function::value;
  virtual Real value(Real t, const Point & p) const override;

};