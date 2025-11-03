//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Material.h"

/**
 * Computes a (possibly damage-dependent) Biot coefficient
 *   alpha = 1 - K(d) / K_s
 * where the damaged bulk modulus K(d) is obtained from the intact modulus K_0
 * through the degradation function g(d) = (1 - d)^2, with d the phase-field
 * (damage) variable. K_s is the solid grain bulk modulus provided as input.
 *
 * Provides MaterialProperty<Real>s named "biot_coefficient_damaged"
 */
class ElkPorousFlowDamagedBiotCoefficient : public Material
{
public:
  static InputParameters validParams();

  ElkPorousFlowDamagedBiotCoefficient(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  /// Damage / phase-field variable
  const VariableValue & _damage;

  /// Intact solid bulk modulus K_0 (computed from the supplied compliance)
  const Real _bulk_modulus_intact;

  /// Solid grain bulk modulus K_s
  const Real _grain_bulk_modulus;

  /// Minimum allowable value for the degradation function g(c)
  const Real _min_degradation;

  /// Damaged Biot coefficient alpha(c)
  MaterialProperty<Real> & _biot_coefficient_damaged;
};
