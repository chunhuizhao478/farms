//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "PorousFlowMaterialVectorBase.h"

/**
 * Computes the Biot modulus using the damage-dependent relationships:
 *   alpha(c) = 1 - K_c / K_s ,
 *   phi(c)   = phi_0 + (1 - phi_0)[1 - (1 - c)^2],
 *   1/M(c)   = phi(c) / K_f + (alpha(c) - phi(c)) / K_s ,
 * where c is the phase-field (damage) variable, K_c is the degraded drained
 * bulk modulus, K_s is the solid grain bulk modulus, and K_f is the fluid bulk modulus.
 */
class ElkPorousFlowDamagedBiotModulus : public PorousFlowMaterialVectorBase
{
public:
  static InputParameters validParams();

  ElkPorousFlowDamagedBiotModulus(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties() override;
  virtual void computeQpProperties() override;

  /// Base-value Biot coefficient (used when damaged coefficient is not supplied)
  const Real _biot_coefficient_const;

  /// Toggle to read biot from material property instead of constant
  const bool _use_damaged_biot;

  /// Damaged Biot coefficient material property (if enabled)
  const MaterialProperty<Real> * _biot_coefficient_damaged_matprop;

  /// Fluid bulk modulus
  const Real _fluid_bulk_modulus;

  /// Solid grain bulk modulus (K_s in the reference model)
  const Real _grain_bulk_modulus;

  /// Toggle to read porosity from material property instead of constant
  const bool _use_damaged_porosity;

  /// Damaged porosity material property (if enabled)
  const MaterialProperty<Real> * _porosity_damaged_matprop;

  /// Constant porosity value (used when use_damaged_porosity=false)
  const Real _porosity_const;

  /// Computed Biot modulus
  MaterialProperty<Real> & _biot_modulus;

  /// Small floor for the denominator when inverting the storativity relation
  const Real _denominator_floor = 1e-20;
};
