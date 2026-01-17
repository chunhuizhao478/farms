//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ADMaterial.h"

/**
 * AD version of ElkPorousFlowDamagedBiotModulus.
 * Computes the Biot modulus using the damage-dependent relationships:
 *   1/M(d) = phi(d) / K_f + (alpha(d) - phi(d)) / K_s
 * where d is the phase-field (damage) variable, alpha is the damaged Biot coefficient,
 * phi is the damaged porosity, K_s is the solid grain bulk modulus, and K_f is the fluid bulk modulus.
 */
class ElkADPorousFlowDamagedBiotModulus : public ADMaterial
{
public:
  static InputParameters validParams();

  ElkADPorousFlowDamagedBiotModulus(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  /// Toggle to read biot from material property instead of constant
  const bool _use_damaged_biot;

  /// Base-value Biot coefficient (used when damaged coefficient is not supplied)
  const Real _biot_coefficient_const;

  /// Damaged Biot coefficient material property (if enabled)
  const ADMaterialProperty<Real> * _biot_coefficient_damaged_matprop;

  /// Fluid bulk modulus
  const Real _fluid_bulk_modulus;

  /// Solid grain bulk modulus (K_s in the reference model)
  const Real _grain_bulk_modulus;

  /// Toggle to read porosity from material property instead of constant
  const bool _use_damaged_porosity;

  /// Damaged porosity material property (if enabled)
  const ADMaterialProperty<Real> * _porosity_damaged_matprop;

  /// Constant porosity value (used when use_damaged_porosity=false)
  const Real _porosity_const;

  /// Computed Biot modulus (AD property for automatic differentiation)
  ADMaterialProperty<Real> & _biot_modulus;

  /// Small floor for the denominator when inverting the storativity relation
  const Real _denominator_floor = 1e-20;
};
