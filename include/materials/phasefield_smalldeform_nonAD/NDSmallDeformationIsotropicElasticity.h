//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#pragma once

#include "NDSmallDeformationElasticityModel.h"
#include "DerivativeMaterialPropertyNameInterface.h"

class NDSmallDeformationIsotropicElasticity : public NDSmallDeformationElasticityModel,
                                            public DerivativeMaterialPropertyNameInterface
{
public:
  static InputParameters validParams();

  NDSmallDeformationIsotropicElasticity(const InputParameters & parameters);

  virtual RankTwoTensor computeStress(const RankTwoTensor & strain) override;

  // Compute Jacobian of the stress w/r/t strain
  virtual RankFourTensor computeJacobian(const RankTwoTensor & strain) override;

protected:
private:
  // @{ Decomposition methods
  virtual RankTwoTensor computeStressNoDecomposition(const RankTwoTensor & strain);
  virtual RankTwoTensor computeStressSpectralDecomposition(const RankTwoTensor & strain);
  virtual RankTwoTensor computeStressVolDevDecomposition(const RankTwoTensor & strain);
  // @}

  //add jacobain of the stress w/r/t strain with decomposition methods
  virtual RankFourTensor computeJacobianNoDecomposition(const RankTwoTensor & strain);
  virtual RankFourTensor computeJacobianSpectralDecomposition(const RankTwoTensor & strain);
  virtual RankFourTensor computeJacobianVolDevDecomposition(const RankTwoTensor & strain);

  // @{ add additional functions for porous flow coupling
  virtual void computeCrackStrainAndOrientation(RealVectorValue & strain_in_crack_dir);
  virtual void updatePermeabilityForCracking();
  // @}

  // @{ Helper functions
  Real Macaulay(const Real x, const bool deriv = false);
  std::vector<Real> Macaulay(const std::vector<Real> & v, const bool deriv = false);
  RankTwoTensor spectralDecomposition(const RankTwoTensor & r2t);
  // @}

  // Compute g and its derivatives
  void computeGDerivatives();

  /// The bulk modulus
  const MaterialProperty<Real> & _K;

  /// The shear modulus
  const MaterialProperty<Real> & _G;

  /// Name of the phase-field variable
  const VariableValue & _d;

  // @{ Strain energy density and its derivative w/r/t damage
  MaterialProperty<Real> & _psie;
  MaterialProperty<Real> & _psie_active;
  MaterialProperty<Real> & _dpsie_dd;
  // @}

  // @{ The degradation function and its derivative w/r/t damage
  MaterialProperty<Real> & _g;
  MaterialProperty<Real> & _dg_dd;
  MaterialProperty<Real> & _d2g_dd2;
  // @}

  // Model type
  const std::string _model_type;

  // Constants
  const Real _eta;

  /// Decomposittion types
  const enum class Decomposition { none, spectral, voldev } _decomposition;

  /// Add additional material properties for porous flow coupling
  //@{ Rotation tensor used to rotate tensors into crack local coordinates
  MaterialProperty<RankTwoTensor> & _crack_rotation;
  const MaterialProperty<RankTwoTensor> & _crack_rotation_old;
  ///@}

  /// @brief define the effective permeability
  MaterialProperty<RealTensorValue> & _effective_perm;
  const MaterialProperty<RealTensorValue> & _effective_perm_old;  

  const bool _porous_flow_coupling; // flag to indicate if porous flow coupling is enabled
  const Real _intrinsic_permeability;
  const Real _coeff_b; // coefficient for the exponential function in the effective permeability

};
