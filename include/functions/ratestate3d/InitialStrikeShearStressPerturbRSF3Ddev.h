/*
Define Function for Initial Shear Stress along Strike Direction
*/

#pragma once

#include "Function.h"

class InitialStrikeShearStressPerturbRSF3Ddev : public Function
{
public:
  InitialStrikeShearStressPerturbRSF3Ddev(const InputParameters & parameters);

  static InputParameters validParams();

  using Function::value;
  virtual Real value(Real t, const Point & p) const override;

  Real _stress_perturb;

};