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
#include "SmearedCrackSofteningBase.h"
#include "Function.h"

/*
AD Farms Compute Smeared Cracking Stress Model
Created by Chunhui Zhao, Apr 20th, 2025
Rewrite the smeared crack model, add energy regularization
Regularization takes place on equivalent strain
Uses automatic differentiation

- Pure Solid Mechanics
- Take regularizated equivalent strain as input

*/

/**
 * ADFarmsComputeSmearedCrackingStressGrads computes the stress for a finite strain
 * material with smeared cracking using automatic differentiation
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
  void computeCrackStrainAndOrientation(RealVectorValue & strain_in_crack_dir);

  // @{ add additional functions for porous flow coupling
  virtual void updatePermeabilityForCracking();
  // @}

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
  const Real _paramA;
  const Real _paramB;
  ///@}

  ///initial damage for crack_damage material property
  const ADVariableValue & _initial_crack_damage;

  //porous flow coupling related parameters
  const bool _porous_flow_coupling; // flag to indicate if porous flow coupling is enabled
  const Real _intrinsic_permeability;

  /// @brief define the effective permeability
  ADMaterialProperty<RealTensorValue> & _effective_perm;
  const MaterialProperty<RealTensorValue> & _effective_perm_old; 

  // Exponential permeability model
  const bool _exponential_permeability_model; // flag to indicate if exponential permeability model is used
  const Real _coeff_b; // coefficient for the exponential function in the effective permeability

  // Darcy-Poiseuille permeability model
  const bool _darcy_poiseuille_permeability_model; // flag to indicate if Darcy-Poiseuille permeability model is used
  const Real _wc; // characteristic width for the Darcy-Poiseuille model
  const Real _perm_exponent; // exponent for the Darcy-Poiseuille model
};