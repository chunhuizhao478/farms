/*
Define Function for Initial Shear Stress Hetergeneoity in XY Direction 
Gaussian Distribution
Chunhui Zhao
*/

#pragma once

#include "Function.h"

class InitialStressXYcontmfbfs3D : public Function
{
public:
  InitialStressXYcontmfbfs3D(const InputParameters & parameters);

  static InputParameters validParams();

  using Function::value;
  virtual Real value(Real t, const Point & p) const override;

  Real _maximum_value;
  Real _length_z;
  Real _nucl_loc_x;
  Real _nucl_loc_y;
  Real _nucl_loc_z;
  Real _nucl_patch_size;

};