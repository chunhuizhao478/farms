//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ComputeGeneralStressBase.h"
#include "SymmetricRankFourTensor.h"
#include "RankTwoTensor.h"
#include "Function.h"

/**
 * Small-strain variant of the spectral smeared-crack model with energy bookkeeping
 * and optional gradient-damage coupling on the equivalent strain.
 *
 * Strain is sourced from the ComputeSmallStrain pipeline via _mechanical_strain.
 */
class FarmsComputeSmearedCrackingStressGradsSpectralSmallStrain : public ComputeGeneralStressBase
{
public:
  static InputParameters validParams();

  FarmsComputeSmearedCrackingStressGradsSpectralSmallStrain(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties() override;
  virtual void computeQpStress() override;

  // Utilities
  Real Macaulay(const Real x, const bool deriv);
  RankTwoTensor spectralPositivePart(const RankTwoTensor & A) const;
  Real equivalentStrainFromTensor(const RankTwoTensor & eps) const;
  RankTwoTensor symmetricBasis(unsigned int i, unsigned int j) const;
  void equivalentStrainDerivativesFD(const RankTwoTensor & eps,
                                     Real delta,
                                     RankTwoTensor & grad,
                                     RankFourTensor & hess) const;

  // Compute principal strains and crack orientation (stores in _crack_rotation)
  void computeCrackStrainAndOrientation(RealVectorValue & strain_in_crack_dir);

  // HM coupling helpers
  virtual void updateSolidBulkCompliance();
  virtual void updatePermeabilityForCracking();

protected:
  // Elasticity tensor (intact)
  const MaterialProperty<RankFourTensor> & _elasticity_tensor;

  // old mechanical strain for proper energy increment dE = E^{n+1} - E^{n}
  const MaterialProperty<RankTwoTensor> & _mechanical_strain_old;

  // Stress old for trapezoidal work
  const MaterialProperty<RankTwoTensor> & _stress_old;

  // Strain energy density outputs
  MaterialProperty<Real> & _psie;
  MaterialProperty<Real> & _psie_active;

  // Energy bookkeeping
  MaterialProperty<Real> & _accumulated_elastic_energy;
  const MaterialProperty<Real> & _accumulated_elastic_energy_old;
  MaterialProperty<Real> & _instant_elastic_energy;
  MaterialProperty<Real> & _fracture_energy;

  // Input and state for damage
  const VariableValue & _cracking_stress;
  MaterialProperty<Real> & _crack_damage;
  const MaterialProperty<Real> & _crack_damage_old;

  // Equivalent strains
  MaterialProperty<Real> & _eqstrain_local;
  const MaterialProperty<Real> & _eqstrain_local_old;
  const VariableValue & _eqstrain_nonlocal;
  MaterialProperty<Real> & _kappa;
  const MaterialProperty<Real> & _kappa_old;

  // Crack orientation
  MaterialProperty<RankTwoTensor> & _crack_rotation;
  const MaterialProperty<RankTwoTensor> & _crack_rotation_old;

  // Damage evolution parameters
  Real _paramA;
  Real _paramB;
  const VariableValue & _initial_crack_damage;

  // Permeability coupling
  const bool _porous_flow_coupling;
  const Real _intrinsic_permeability;
  MaterialProperty<RealTensorValue> & _effective_perm;
  const MaterialProperty<RealTensorValue> & _effective_perm_old;
  const bool _exponential_permeability_model;
  const Real _coeff_b;
  const bool _darcy_poiseuille_permeability_model;
  const Real _wc;
  const Real _perm_exponent;

  // Damaged solid bulk compliance C_s(d)
  MaterialProperty<Real> & _solid_bulk_compliance_damaged;

  // Gradient-damage parameters
  Real _h;
  Real _fd_delta;
};
