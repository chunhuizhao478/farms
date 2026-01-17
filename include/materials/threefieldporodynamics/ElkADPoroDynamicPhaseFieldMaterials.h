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
 * ADMaterial used in three-field poro-dynamics phase-field simulations.
 * This material assembles all damage-dependent hydraulic properties:
 * - biot_coefficient (from ElkADPorousFlowDamagedBiotCoefficient)
 * - porosity (from ElkADPorousFlowDamagedPorosity)
 * - biot_modulus (from ElkADPorousFlowDamagedBiotModulus)
 * - effective_perm (from phase-field elasticity model)
 *
 * And computes:
 * - density = rho_s * (1 - phi) + rho_f * phi
 * - volumetric strain from elastic strain tensor
 */
class ElkADPoroDynamicPhaseFieldMaterials : public ADMaterial
{
public:
  static InputParameters validParams();

  ElkADPoroDynamicPhaseFieldMaterials(const InputParameters & parameters);

  virtual void initQpStatefulProperties() override;
  virtual void computeQpProperties() override;

protected:

  /// Solid density value
  const Real _rhos_val;

  /// Fluid density value
  const Real _rhof_val;

  /// Tortosity value
  const Real _tortosity_val;

  /// Viscosity value
  const Real _viscosity_val;

  /// Solid grain bulk modulus (for biot modulus calculation if not using material property)
  const Real _grain_bulk_modulus;

  /// Fluid bulk modulus (for biot modulus calculation if not using material property)
  const Real _fluid_bulk_modulus;

  /// Toggle to use damaged properties from material
  const bool _use_damaged_properties;

  /// Fallback constant values if not using damaged properties
  const Real _porosity_const;
  const Real _biot_coefficient_const;

  /// Damaged Biot coefficient material property (if enabled)
  const ADMaterialProperty<Real> * _biot_coefficient_damaged;

  /// Damaged porosity material property (if enabled)
  const ADMaterialProperty<Real> * _porosity_damaged;

  /// Damaged Biot modulus material property (if enabled)
  const ADMaterialProperty<Real> * _biot_modulus_damaged;

  /// Effective permeability tensor (from phase-field elasticity model)
  const ADMaterialProperty<RankTwoTensor> & _effective_perm;

  /// Solid elastic strain (for volumetric strain computation)
  const ADMaterialProperty<RankTwoTensor> & _elastic_strain;

  // Output material properties
  /// Material property fluid density
  ADMaterialProperty<Real> & _rhof;

  /// Material property density
  ADMaterialProperty<Real> & _rho;

  /// Material property porosity
  ADMaterialProperty<Real> & _porosity;

  /// Material property tortosity
  ADMaterialProperty<Real> & _tortosity;

  /// Material property viscosity
  ADMaterialProperty<Real> & _viscosity;

  /// Material property biot modulus
  ADMaterialProperty<Real> & _biot_modulus;

  /// Material property biot coefficients
  ADMaterialProperty<Real> & _biot_coefficient;

  /// Material property permeability (tensor)
  ADMaterialProperty<RankTwoTensor> & _permeability;

  /// Material property volumetric strain
  ADMaterialProperty<Real> & _vol_strain;
};
