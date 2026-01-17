//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "SmallDeformationElasticityModel.h"
#include "DerivativeMaterialPropertyNameInterface.h"

/**
 * AD version of isotropic elasticity with phase-field damage coupling.
 * Provides:
 * - Stress computation with spectral/volumetric-deviatoric decomposition
 * - Degradation function g(d) and derivatives (computed internally)
 * - Strain energy density and derivatives for phase-field driving force
 * - Damage-dependent permeability (Darcy-Poiseuille model)
 *
 * This follows the non-AD NDSmallDeformationIsotropicElasticity but uses AD types.
 * Compatible with three-field poroelastodynamics formulation.
 */
class ADSmallDeformationIsotropicElasticityPF : public SmallDeformationElasticityModel,
                                                 public DerivativeMaterialPropertyNameInterface
{
public:
  static InputParameters validParams();

  ADSmallDeformationIsotropicElasticityPF(const InputParameters & parameters);

  virtual ADRankTwoTensor computeStress(const ADRankTwoTensor & strain) override;

protected:
  virtual void initQpStatefulProperties() override;

  /// Compute stress based on decomposition type
  ADRankTwoTensor computeStressNoDecomposition(const ADRankTwoTensor & strain);
  ADRankTwoTensor computeStressSpectralDecomposition(const ADRankTwoTensor & strain);
  ADRankTwoTensor computeStressVolDevDecomposition(const ADRankTwoTensor & strain);

  /// Compute degradation function and derivatives
  void computeGDerivatives();

  /// Update permeability based on damage (for porous flow coupling)
  void updatePermeabilityForCracking();

  /// Compute crack strain and orientation
  void computeCrackStrainAndOrientation(RealVectorValue & strain_in_crack_dir);

  /// Helper: Macaulay bracket (AD version)
  ADReal Macaulay(const ADReal & x);
  std::vector<ADReal> Macaulay(const std::vector<ADReal> & v);

  /// Spectral decomposition helper (AD version)
  ADRankTwoTensor spectralDecomposition(const ADRankTwoTensor & r2t);

  /// The bulk modulus material property
  const ADMaterialProperty<Real> & _K;

  /// The shear modulus material property
  const ADMaterialProperty<Real> & _G;

  /// The phase-field (damage) variable
  const ADVariableValue & _d;

  /// Model type (AT1, AT2)
  const std::string _model_type;

  /// Minimum degradation (residual stiffness)
  const Real _eta;

  /// Decomposition type
  const enum class Decomposition { none, spectral, voldev } _decomposition;

  /// Porous flow coupling enabled
  const bool _porous_flow_coupling;

  /// Intrinsic permeability
  const Real _intrinsic_permeability;

  /// Darcy-Poiseuille permeability model parameters
  const bool _darcy_poiseuille_permeability_model;
  const Real _wc;
  const Real _perm_exponent;

  /// Exponential permeability model parameters
  const bool _exponential_permeability_model;
  const Real _coeff_b;

  // Output material properties

  /// Strain energy density
  ADMaterialProperty<Real> & _psie;

  /// Active (tensile) strain energy density
  ADMaterialProperty<Real> & _psie_active;

  /// Inactive (compressive) strain energy density
  ADMaterialProperty<Real> & _psie_inactive;

  /// Derivative of strain energy density w.r.t. damage
  ADMaterialProperty<Real> & _dpsie_dd;

  /// Degradation function g(d)
  ADMaterialProperty<Real> & _g;

  /// First derivative of degradation function
  ADMaterialProperty<Real> & _dg_dd;

  /// Second derivative of degradation function
  ADMaterialProperty<Real> & _d2g_dd2;

  /// Crack rotation tensor (for permeability)
  ADMaterialProperty<RankTwoTensor> & _crack_rotation;
  const MaterialProperty<RankTwoTensor> & _crack_rotation_old;

  /// Effective permeability (RankTwoTensor for AD compatibility with three-field)
  ADMaterialProperty<RankTwoTensor> & _effective_perm;
  const MaterialProperty<RankTwoTensor> & _effective_perm_old;
};
