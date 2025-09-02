//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ColumnMajorMatrix.h"
#include "ComputeMultipleInelasticStress.h"
#include "SmearedCrackSofteningBase.h"
#include "Function.h"

/*
Farms Compute Smeared Cracking Stress Model
Created by Chunhui Zhao, Apr 20th, 2025
Rewrite the smeared crack model, add energy regularization
Regularization takes place on equivalent strain

- Pure Solid Mechanics
- Take regularizated equivalent strain as input

*/

/**
 * FarmsComputeSmearedCrackingStressGradsSpectral computes the stress for a finite strain
 * material with smeared cracking
 */
class FarmsComputeSmearedCrackingStressGradsSpectral : public ComputeMultipleInelasticStress
{
public:
  static InputParameters validParams();

  FarmsComputeSmearedCrackingStressGradsSpectral(const InputParameters & parameters);

  virtual void initQpStatefulProperties() override;
  virtual void computeQpStress() override;

protected:

  /**
   * Compute the crack strain in the crack coordinate system. Also
   * computes the crack orientations, and stores in _crack_rotation.
   * @param strain_in_crack_dir Computed strains in crack directions
   */
  void computeCrackStrainAndOrientation(RealVectorValue & strain_in_crack_dir);

  // Update solid bulk compliance
  virtual void updateSolidBulkCompliance();

  // @{ add additional functions for porous flow coupling
  virtual void updatePermeabilityForCracking();
  // @}

  // Macaulay bracket function
  Real Macaulay(const Real x, const bool deriv);

  // Spectral decomposition for positive part of a symmetric tensor
  RankTwoTensor spectralPositivePart(const RankTwoTensor & A) const;

  // Equivalent strain from a symmetric tensor (Mazars-type: sqrt(sum <eps_i>^2))
  Real equivalentStrainFromTensor(const RankTwoTensor & eps) const;

  // Build symmetric basis tensor E^(ij) with 1 on diagonal, 1/2 on symmetric off-diagonals
  RankTwoTensor symmetricBasis(unsigned int i, unsigned int j) const;

  // Finite-difference derivatives of equivalent strain wrt strain tensor
  // Returns gradient G such that de ≈ G : dε, and Hessian H such that
  // d^2 e ≈ (H :: (dε ⊗ dε)) where :: is double contraction on both pairs
  void equivalentStrainDerivativesFD(const RankTwoTensor & eps,
                                     Real delta,
                                     RankTwoTensor & grad,
                                     RankFourTensor & hess) const;

  // Analytical derivatives using spectral projectors (preferred)
  void equivalentStrainDerivativesAnalytical(const RankTwoTensor & eps,
                                             RankTwoTensor & grad,
                                             RankFourTensor & hess) const;

  // @{ Strain energy density outputs (degraded)
  // psie_active is the tensile (active) part of the intact energy
  // psie = g * psie_active + psie_inactive, where g = 1 - omega
  MaterialProperty<Real> & _psie;
  MaterialProperty<Real> & _psie_active;
  ///@}

  // @{ Elastic energy bookkeeping for dissipation accounting
  // Accumulated elastic energy Ea (history) and its old value
  MaterialProperty<Real> & _accumulated_elastic_energy;
  const MaterialProperty<Real> & _accumulated_elastic_energy_old;
  // Instant elastic energy Ei = 1/2 * sigma : epsilon
  MaterialProperty<Real> & _instant_elastic_energy;
  // Fracture (dissipated) energy = Ea - Ei
  MaterialProperty<Real> & _fracture_energy;
  ///@}

  ///@{ Input parameters for smeared crack models

  /// Threshold at which cracking initiates if tensile stress exceeds it
  const VariableValue & _cracking_stress;

  //@{ Damage (goes from 0 to 1) in crack directions
  //Damage is treated as a scalar variable
  MaterialProperty<Real> & _crack_damage;
  const MaterialProperty<Real> & _crack_damage_old;
  ///@}

  ///@{equivalent local strain value
  MaterialProperty<Real> & _eqstrain_local;
  const MaterialProperty<Real> & _eqstrain_local_old;
  ///@}

  ///@{equivalent nonlocal strain value
  const VariableValue & _eqstrain_nonlocal;
  ///@}

  ///@{maximum equivalent strain value
  MaterialProperty<Real> & _kappa;
  const MaterialProperty<Real> & _kappa_old;
  ///@}

  //@{ Rotation tensor used to rotate tensors into crack local coordinates
  MaterialProperty<RankTwoTensor> & _crack_rotation;
  const MaterialProperty<RankTwoTensor> & _crack_rotation_old;
  ///@}

  ///@{ Parameters for the damage evolution law
  Real _paramA;
  Real _paramB;
  ///@}

  ///initial damage for crack_damage material property
  const VariableValue & _initial_crack_damage;

  /// Vector helper to update local elasticity tensor
  std::vector<Real> _local_elastic_vector;
  /// Variables used by multiple methods within the calculation for a single material point
  RankFourTensor _local_elasticity_tensor;

  //porous flow coupling related parameters
  const bool _porous_flow_coupling; // flag to indicate if porous flow coupling is enabled
  const Real _intrinsic_permeability;

  /// @brief define the effective permeability
  MaterialProperty<RealTensorValue> & _effective_perm;
  const MaterialProperty<RealTensorValue> & _effective_perm_old; 

  // Exponential permeability model
  const bool _exponential_permeability_model; // flag to indicate if exponential permeability model is used
  const Real _coeff_b; // coefficient for the exponential function in the effective permeability

  // Darcy-Poiseuille permeability model
  const bool _darcy_poiseuille_permeability_model; // flag to indicate if Darcy-Poiseuille permeability model is used
  const Real _wc; // characteristic width for the Darcy-Poiseuille model
  const Real _perm_exponent; // exponent for the Darcy-Poiseuille model

  // Damaged solid bulk compliance C_s(d) = 1 / (g(d) * K)
  MaterialProperty<Real> & _solid_bulk_compliance_damaged;

  // Gradient damage coupling modulus (energy term 1/2 h (e - e~)^2)
  Real _h;
  // Finite-difference step for equivalent strain derivatives
  Real _fd_delta;

};