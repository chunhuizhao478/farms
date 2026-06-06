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
#include "RankTwoTensor.h"

/**
 * Computes the damaged porosity stored in the property
 *   PorousFlow_porosity_{qp,nodal}_damaged
 * using one of two selectable laws (porosity_update_model):
 *
 *   - "damage" (default): damage-driven, bounded-maximum porosity
 *       phi(d) = phi_0 + (1 - phi_0)[1 - (1 - d)^2]
 *     which saturates at phi = 1 and is capped by porosity_upper_bound.
 *
 *   - "strain": strain-based update of Liu et al. (2024, CMAME 429:117165) eq. (40)
 *       phi(eps) = phi_0 + eps_1
 *     where eps_1 is the maximum (most-tensile) principal value of the kinematic
 *     strain tensor (strain_property, default "mechanical_strain"). This branch is
 *     only available at quadrature points and is evaluated instantaneously (no
 *     history): if the strain relaxes the porosity decreases again.
 *
 * Both laws are clamped to [porosity_lower_bound, porosity_upper_bound]. The
 * undamaged PorousFlow porosity material (e.g., PorousFlowPorosityConst) should
 * still be supplied separately when its values are required.
 */
class ElkPorousFlowDamagedPorosity : public Material
{
public:
  static InputParameters validParams();

  ElkPorousFlowDamagedPorosity(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  /// Which porosity-update law to evaluate
  enum class PorosityUpdateModel
  {
    DAMAGE,
    STRAIN
  };
  const PorosityUpdateModel _porosity_update_model;

  /// Damage / phase-field variable (used by the DAMAGE model; unused by STRAIN)
  const VariableValue & _damage;

  /// Initial porosity phi_0
  const Real _initial_porosity;

  /// Minimum and maximum porosity bounds applied after the update
  const Real _porosity_lower_bound;
  const Real _porosity_upper_bound;

  /// Kinematic strain tensor; bound only when _porosity_update_model == STRAIN
  /// (nullptr otherwise). eps_1 = max principal value drives the eq. (40) update.
  const MaterialProperty<RankTwoTensor> * _mechanical_strain;

  /// Damaged porosity property
  MaterialProperty<Real> & _porosity_damaged;
};
