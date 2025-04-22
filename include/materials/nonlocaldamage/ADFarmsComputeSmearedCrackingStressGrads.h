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
#include "ADComputeMultipleInelasticStress.h"
#include "ADSmearedCrackSofteningBase.h"
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
 * ADFarmsComputeSmearedCrackingStressGrads computes the stress for a finite strain
 * material with smeared cracking
 */
class ADFarmsComputeSmearedCrackingStressGrads : public ADComputeMultipleInelasticStress
{
public:
  static InputParameters validParams();

  ADFarmsComputeSmearedCrackingStressGrads(const InputParameters & parameters);

  virtual void initQpStatefulProperties() override;
  virtual void computeQpStress() override;

protected:

  /**
   * Compute the crack strain in the crack coordinate system. Also
   * computes the crack orientations, and stores in _crack_rotation.
   * @param strain_in_crack_dir Computed strains in crack directions
   */
  void computeCrackStrainAndOrientation(ADRealVectorValue & strain_in_crack_dir);

  /**
   * Update the local elasticity tensor (_local_elasticity_tensor)
   * due to the effects of cracking.
   */
  // void updateLocalElasticityTensor();

  ///@{ Input parameters for smeared crack models

  /// Threshold at which cracking initiates if tensile stress exceeds it
  const ADVariableValue & _cracking_stress;

  //@{ Damage (goes from 0 to 1) in crack directions
  //Damage is treated as a scalar variable
  ADMaterialProperty<Real> & _crack_damage;
  const MaterialProperty<Real> & _crack_damage_old;
  ///@}

  ///@{equivalent local strain value
  ADMaterialProperty<Real> & _eqstrain_local;
  const MaterialProperty<Real> & _eqstrain_local_old;
  ///@}

  ///@{equivalent nonlocal strain value
  const ADVariableValue & _eqstrain_nonlocal;
  ///@}

  ///@{maximum equivalent strain value
  ADMaterialProperty<Real> & _kappa;
  const MaterialProperty<Real> & _kappa_old;
  ///@}

  //@{ Rotation tensor used to rotate tensors into crack local coordinates
  ADMaterialProperty<RankTwoTensor> & _crack_rotation;
  const MaterialProperty<RankTwoTensor> & _crack_rotation_old;
  ///@}

  ///@{ Parameters for the damage evolution law
  ADReal _paramA;
  ADReal _paramB; 
  ///@}

};