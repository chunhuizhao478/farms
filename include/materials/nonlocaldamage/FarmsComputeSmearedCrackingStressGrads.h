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
 * FarmsComputeSmearedCrackingStressGrads computes the stress for a finite strain
 * material with smeared cracking
 */
class FarmsComputeSmearedCrackingStressGrads : public ComputeMultipleInelasticStress
{
public:
  static InputParameters validParams();

  FarmsComputeSmearedCrackingStressGrads(const InputParameters & parameters);

  virtual void initQpStatefulProperties() override;
  virtual void computeQpStress() override;

protected:

  /**
   * Compute the crack strain in the crack coordinate system. Also
   * computes the crack orientations, and stores in _crack_rotation.
   * @param strain_in_crack_dir Computed strains in crack directions
   */
  void computeCrackStrainAndOrientation(RealVectorValue & strain_in_crack_dir);

  /**
   * Update the local elasticity tensor (_local_elasticity_tensor)
   * due to the effects of cracking.
   */
  // void updateLocalElasticityTensor();

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

};