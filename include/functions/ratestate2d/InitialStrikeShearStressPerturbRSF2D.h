/*
Define Function for Initial Shear Stress along Strike Direction
*/

#pragma once

#include "Function.h"

class InitialStrikeShearStressPerturbRSF2D : public Function
{
public:
  InitialStrikeShearStressPerturbRSF2D(const InputParameters & parameters);

  static InputParameters validParams();

  using Function::value;
  virtual Real value(Real t, const Point & p) const override;

};