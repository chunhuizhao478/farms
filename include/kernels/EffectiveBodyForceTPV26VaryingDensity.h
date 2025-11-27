/*
Kernel to apply effective body force accounting for pore pressure gradient
with depth-varying density from TPV32 profile.
The effective body force is: f_eff = rho(z) * g - dPf/dz
This extends EffectiveBodyForceTPV26 to support depth-varying density.
Created By Chunhui Zhao, Nov 26th, 2025
*/

#pragma once

#include "Kernel.h"
#include <vector>

class EffectiveBodyForceTPV26VaryingDensity : public Kernel
{
public:
  static InputParameters validParams();

  EffectiveBodyForceTPV26VaryingDensity(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;

private:
  // TPV32 piecewise-linear density (kg/m^3) at depth d (m)
  Real tpv32_density_at_depth(Real d) const;

  Real _fluid_density; //fluid density in kg/m^3
  Real _gravity; //gravity in m/s^2
  bool _use_overpressure; //flag to use overpressure in the calculation
  Real _overpressure_depth_A; //depth at which overpressure starts to transition
  Real _overpressure_depth_B; //depth at which overpressure stops to transition
  Real _overpressure_rho_ref; //reference rock density for overpressure (default: TPV32 final value = 2900)
};
