/*
Define Function for Initial Static Friction Coefficient for benchmark
*/

#pragma once

#include "Function.h"

class SolutionUserObjectBase;

class InitialShearStressCDBM : public Function
{
public:
  InitialShearStressCDBM(const InputParameters & parameters);

  static InputParameters validParams();

  using Function::value;
  virtual Real value(Real t, const Point & p) const override;

  Real _peak_value;

  Real _nucl_center_x;
  Real _nucl_center_z;

  Real _nucl_size;
  
  Real _elem_size;

  Real _domain_value;

};