/*
Define Function for Initial Shear Stress Hetergeneoity in XY Direction 
Will be Passed directly to "ComputeEigenstrainFromInitialStress" as variables
Chunhui Zhao

cluster run stepover_3
*/

#pragma once

#include "Function.h"

class InitialStressXYfinitelen_multiactivation : public Function
{
public:
  InitialStressXYfinitelen_multiactivation(const InputParameters & parameters);

  static InputParameters validParams();

  using Function::value;
  virtual Real value(Real t, const Point & p) const override;

};