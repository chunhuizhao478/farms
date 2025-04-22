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
Farms Crack Band Model
Created by Chunhui Zhao, Apr 18th, 2025
Rewrite the smeared crack model, add energy regularization
Regularization takes place on equivalent strain

- Pure Solid Mechanics
- Take regularizated equivalent strain as input

*/

/**
 * FarmsCrackBandModel computes the stress for a finite strain
 * material with smeared cracking
 */
class FarmsCrackBandModel : public ComputeMultipleInelasticStress
{
public:
  static InputParameters validParams();

  FarmsCrackBandModel(const InputParameters & parameters);

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
   * Compute the crack strain in the crack coordinate system. Also
   * computes the crack orientations, and stores in _crack_rotation.
   * @param stress_in_crack_dir Computed strains in crack directions
   */
  void computeCrackStressAndOrientation(RealVectorValue & stress_in_crack_dir);

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

  ///@{equivalent strain value
  MaterialProperty<Real> & _eqstrain;
  const MaterialProperty<Real> & _eqstrain_old;
  ///@}

  ///@{maximum equivalent strain value
  MaterialProperty<Real> & _kappa;
  const MaterialProperty<Real> & _kappa_old;
  ///@}

  //@{ Rotation tensor used to rotate tensors into crack local coordinates
  MaterialProperty<RankTwoTensor> & _crack_rotation;
  const MaterialProperty<RankTwoTensor> & _crack_rotation_old;
  ///@}

  ///@{parameters for crack band model
  Real _hb; //Crack‐band size (element length)
  Real _Gf; //Fracture energy
  ///@}

};