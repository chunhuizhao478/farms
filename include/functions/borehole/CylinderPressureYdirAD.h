/*
Define Function for Initial Cylinder Pressure Along a Specific Direction
Chunhui Zhao, May 19
*/

#pragma once

#include "Function.h"

class CylinderPressureYdirAD : public Function
{
public:
  CylinderPressureYdirAD(const InputParameters & parameters);

  static InputParameters validParams();

  using Function::value;
  virtual Real value(Real t, const Point & p) const override;

  Real _value;

};