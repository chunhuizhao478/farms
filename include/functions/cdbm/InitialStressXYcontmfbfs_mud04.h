/*
Define Function for Initial Shear Stress Hetergeneoity in XY Direction 
Gaussian Distribution
Chunhui Zhao
*/

#pragma once

#include "Function.h"

class InitialStressXYcontmfbfs_mud04 : public Function
{
public:
  InitialStressXYcontmfbfs_mud04(const InputParameters & parameters);

  static InputParameters validParams();

  using Function::value;
  virtual Real value(Real t, const Point & p) const override;

};