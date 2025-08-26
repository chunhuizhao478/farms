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
 *   alpha = 1 - K / K_s
 * where K is the (damaged) drained bulk modulus of the skeleton, obtained as
 *   K = 1 / solid_bulk_compliance_damaged
 * and K_s is the solid grain bulk modulus provided as input.
 *
 * Provides a MaterialProperty<Real> named "biot_coefficient".
 */
class ElkPorousFlowDamagedBiotCoefficient : public Material
{
public:
  static InputParameters validParams();

  ElkPorousFlowDamagedBiotCoefficient(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  /// Solid grain bulk modulus K_s
  const Real _grain_bulk_modulus;

  /// Damaged solid bulk compliance (1 / K)
  const MaterialProperty<Real> & _solid_bulk_compliance_damaged;

  /// Output Biot coefficient alpha
  MaterialProperty<Real> & _biot_coefficient;
};
