/*
Define Function for Initial Shear Stress along Strike Direction
*/

#pragma once

#include "Function.h"

class InitialStrikeShearStressPerturbRSF2Ddev : public Function
{
public:
  InitialStrikeShearStressPerturbRSF2Ddev(const InputParameters & parameters);

  static InputParameters validParams();

  using Function::value;
  virtual Real value(Real t, const Point & p) const override;

  //stress perturbation value
  Real _stress_perturb;

};