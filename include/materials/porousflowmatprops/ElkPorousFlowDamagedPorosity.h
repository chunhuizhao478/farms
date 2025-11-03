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
 * Computes the damage-dependent porosity
 *   phi(c) = phi_0 + (1 - phi_0)[1 - (1 - c)^2]
 * and stores the result in the property
 *   PorousFlow_porosity_{qp,nodal}_damaged.
 * The undamaged PorousFlow porosity material (e.g., PorousFlowPorosityConst)
 * should still be supplied separately when its values are required.
 */
class ElkPorousFlowDamagedPorosity : public Material
{
public:
  static InputParameters validParams();

  ElkPorousFlowDamagedPorosity(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  /// Damage / phase-field variable
  const VariableValue & _damage;

  /// Initial porosity phi_0
  const Real _initial_porosity;

  /// Minimum and maximum porosity bounds applied after the update
  const Real _porosity_lower_bound;
  const Real _porosity_upper_bound;

  /// Damage-dependent porosity property
  MaterialProperty<Real> & _porosity_damaged;
};
