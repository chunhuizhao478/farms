/*
Define Function for Initial Shear Stress Hetergeneoity in XY Direction 
Gaussian Distribution
Chunhui Zhao
*/

#pragma once

#include "Function.h"

class InitialStresscontmfbfs3D : public Function
{
public:
  InitialStresscontmfbfs3D(const InputParameters & parameters);

  static InputParameters validParams();

  using Function::value;
  virtual Real value(Real t, const Point & p) const override;

  Real _maximum_value;
  Real _length_z;

};