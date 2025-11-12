/*
Kernel to apply effective body force accounting for pore pressure gradient.
The effective body force is: f_eff = rho * g - dPf/dz
This replaces the standard BodyForce kernel when using overpressure models.
Created By Chunhui Zhao, Nov 9th, 2025
*/

#pragma once

#include "Kernel.h"

class EffectiveBodyForceTPV26 : public Kernel
{
public:
  static InputParameters validParams();

  EffectiveBodyForceTPV26(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;

private:
  Real _fluid_density; //fluid density in kg/m^3
  Real _rock_density; //rock density in kg/m^3
  Real _gravity; //gravity in m/s^2
  bool _use_overpressure; //flag to use overpressure in the calculation
  Real _overpressure_depth_A; //depth at which overpressure starts to transition
  Real _overpressure_depth_B; //depth at which overpressure stops to transition
  bool _overpressure_loweffective; //flag to use low effective stress overpressure
  Real _lambda_pp; //pore pressure ratio for low effective stress overpressure
};
