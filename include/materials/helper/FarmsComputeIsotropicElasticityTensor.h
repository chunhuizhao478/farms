//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ComputeElasticityTensorBase.h"

/**
 * FarmsComputeIsotropicElasticityTensor defines an isotropic elasticity tensor
 * where the shear modulus is updated using the initial damage and xi_o via
 * a damage–breakage parameter gamma_r(lambda, mu, xi_o).
 *
 * Effective shear modulus used:
 *   mu_eff = mu_0 + xi_o * initial_damage * gamma_r
 *
 * The lambda parameter is used as provided. This class mirrors the interface
 * of ComputeIsotropicElasticityTensor but requires the pair (lambda, shear_modulus)
 * and an additional xi_o parameter, and reads the 'initial_damage' material property.
 */
template <bool is_ad, typename T>
class FarmsComputeIsotropicElasticityTensorTempl : public ComputeElasticityTensorBaseTempl<is_ad, T>
{
public:
  static InputParameters validParams();

  FarmsComputeIsotropicElasticityTensorTempl(const InputParameters & parameters);

  virtual void initialSetup() override { /* nothing precomputed; depends on qp damage */ }

  virtual void residualSetup() override { /* nothing precomputed; depends on qp damage */ }

protected:
  virtual void computeQpElasticityTensor() override;

  // Inputs
  const Real & _lambda_0;
  const Real & _mu_0;
  const Real & _xi_o;

  // Material property: initial damage (scalar)
  const MaterialProperty<Real> & _initial_damage;

  // Optional material-property inputs for lambda and mu
  const MaterialProperty<Real> * _lambda_input_prop;
  const MaterialProperty<Real> * _mu_input_prop;
  bool _use_input_props;

  // Helper: compute gamma_r(lambda_0, mu_0, xi_o)
  Real computeGammaR(Real lambda0, Real mu0, Real xi0) const;

  // Cached constant gamma_r
  Real _gamma_r;

  using ComputeElasticityTensorBaseTempl<is_ad, T>::_elasticity_tensor;
  using ComputeElasticityTensorBaseTempl<is_ad, T>::_effective_stiffness;
  using ComputeElasticityTensorBaseTempl<is_ad, T>::_qp;
};

typedef FarmsComputeIsotropicElasticityTensorTempl<false, RankFourTensor>
    FarmsComputeIsotropicElasticityTensor;
typedef FarmsComputeIsotropicElasticityTensorTempl<true, RankFourTensor>
    ADFarmsComputeIsotropicElasticityTensor;
typedef FarmsComputeIsotropicElasticityTensorTempl<false, SymmetricRankFourTensor>
    FarmsSymmetricIsotropicElasticityTensor;
typedef FarmsComputeIsotropicElasticityTensorTempl<true, SymmetricRankFourTensor>
    ADFarmsSymmetricIsotropicElasticityTensor;
