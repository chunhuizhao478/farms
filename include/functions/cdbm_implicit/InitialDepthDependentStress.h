/*
Define Function for Initial Depth Dependent Stress
Chunhui Zhao, Mar 19, 2025
*/

#pragma once

#include "Function.h"

class InitialDepthDependentStress : public Function
{
public:
  InitialDepthDependentStress(const InputParameters & parameters);

  static InputParameters validParams();

  using Function::value;
  virtual Real value(Real t, const Point & p) const override;

  Real _i; //index
  Real _j; //index
  bool _sign; //sign of the stress
  Real _fluid_density; //kg/m^3 fluid density
  Real _rock_density; //kg/m^3 rock density
  Real _gravity; //m/s^2
  Real _bxx; //coefficient for sigma_xx
  Real _byy; //coefficient for sigma_yy
  Real _bxy; //coefficient for sigma_xy
  Real _linear_variation_cutoff_distance; //linear variation cutoff distance  

};