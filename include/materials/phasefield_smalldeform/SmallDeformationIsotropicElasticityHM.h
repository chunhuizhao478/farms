//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#pragma once

#include "SmallDeformationElasticityModel.h"
#include "DerivativeMaterialPropertyNameInterface.h"

class SmallDeformationIsotropicElasticityHM : public SmallDeformationElasticityModel,
                                            public DerivativeMaterialPropertyNameInterface
{
public:
  static InputParameters validParams();

  SmallDeformationIsotropicElasticityHM(const InputParameters & parameters);

  virtual ADRankTwoTensor computeStress(const ADRankTwoTensor & strain) override;

protected:
private:
  // @{ Decomposition methods
  virtual ADRankTwoTensor computeStressNoDecomposition(const ADRankTwoTensor & strain);
  virtual ADRankTwoTensor computeStressSpectralDecomposition(const ADRankTwoTensor & strain);
  virtual ADRankTwoTensor computeStressVolDevDecomposition(const ADRankTwoTensor & strain);
  // @}

  // @{ add additional functions for porous flow coupling
  virtual void computeCrackStrainAndOrientation(ADRealVectorValue & strain_in_crack_dir);
  virtual void updatePermeabilityForCracking();
  // @}

  /// The bulk modulus
  const ADMaterialProperty<Real> & _K;

  /// The shear modulus
  const ADMaterialProperty<Real> & _G;

  /// Name of the phase-field variable
  const VariableName _d_name;

  // @{ Strain energy density and its derivative w/r/t damage
  const MaterialPropertyName _psie_name;
  ADMaterialProperty<Real> & _psie;
  ADMaterialProperty<Real> & _psie_active;
  ADMaterialProperty<Real> & _dpsie_dd;
  // @}

  // @{ The degradation function and its derivative w/r/t damage
  const MaterialPropertyName _g_name;
  const ADMaterialProperty<Real> & _g;
  const ADMaterialProperty<Real> & _dg_dd;
  // @}

  /// Decomposittion types
  const enum class Decomposition { none, spectral, voldev } _decomposition;

  /// Add additional material properties for porous flow coupling
  //@{ Rotation tensor used to rotate tensors into crack local coordinates
  ADMaterialProperty<RankTwoTensor> & _crack_rotation;
  const MaterialProperty<RankTwoTensor> & _crack_rotation_old;
  ///@}

  /// @brief define the effective permeability
  ADMaterialProperty<RealTensorValue> & _effective_perm;
  const MaterialProperty<RealTensorValue> & _effective_perm_old;  

  const bool _porous_flow_coupling; // flag to indicate if porous flow coupling is enabled
  const Real _intrinsic_permeability;

  // Exponential permeability model
  const bool _exponential_permeability_model; // flag to indicate if exponential permeability model is used
  const Real _coeff_b; // coefficient for the exponential function in the effective permeability
  // Darcy-Poiseuille permeability model
  const bool _darcy_poiseuille_permeability_model; // flag to indicate if Darcy-Poiseuille permeability model is used
  const Real _wc; // characteristic width for the Darcy-Poiseuille model
  const Real _perm_exponent; // exponent for the Darcy-Poiseuille model

  /// Name of the phase-field variable
  const ADVariableValue & _d;
};
