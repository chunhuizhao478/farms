//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Material.h"

/**
 * RadiationDampingMaterial computes the radiation damping coefficient
 * for quasi-dynamic SEAS simulations.
 *
 * The radiation damping term provides an approximation to the dynamic
 * stress transfer, allowing use of quasi-static elasticity solvers.
 *
 * Radiation damping coefficient: η = μ / (2 * cs)
 *
 * where:
 *   μ = shear modulus
 *   cs = sqrt(μ/ρ) = shear wave speed
 *   ρ = density
 *
 * For BP2 benchmark:
 *   μ = 32.04 GPa
 *   ρ = 2670 kg/m³
 *   cs = 3464 m/s
 *   η = 4.625e6 Pa·s/m
 */
class RadiationDampingMaterial : public Material
{
public:
  static InputParameters validParams();

  RadiationDampingMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  /// Shear modulus
  const Real _shear_modulus;

  /// Density
  const Real _density;

  /// Output: radiation damping coefficient
  MaterialProperty<Real> & _radiation_damping;

  /// Output: shear wave speed
  MaterialProperty<Real> & _shear_wave_speed;
};
