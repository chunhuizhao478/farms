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
 * AD version of ElkPorousFlowDamagedPorosity.
 * Computes the damage-dependent porosity:
 *   phi(d) = phi_0 + (1 - phi_0)[1 - (1 - d)^2]
 * and stores the result in the property:
 *   PorousFlow_porosity_qp_damaged
 */
class ElkADPorousFlowDamagedPorosity : public ADMaterial
{
public:
  static InputParameters validParams();

  ElkADPorousFlowDamagedPorosity(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  /// Damage / phase-field variable
  const ADVariableValue & _damage;

  /// Initial porosity phi_0
  const Real _initial_porosity;

  /// Minimum and maximum porosity bounds applied after the update
  const Real _porosity_lower_bound;
  const Real _porosity_upper_bound;

  /// Damage-dependent porosity property
  ADMaterialProperty<Real> & _porosity_damaged;
};
