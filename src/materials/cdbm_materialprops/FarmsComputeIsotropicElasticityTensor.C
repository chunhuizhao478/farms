//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FarmsComputeIsotropicElasticityTensor.h"
#include "RankFourTensor.h"
#include <cmath>

registerMooseObject("farmsApp", FarmsComputeIsotropicElasticityTensor);
registerMooseObject("farmsApp", ADFarmsComputeIsotropicElasticityTensor);
registerMooseObject("farmsApp", FarmsSymmetricIsotropicElasticityTensor);
registerMooseObject("farmsApp", ADFarmsSymmetricIsotropicElasticityTensor);

template <bool is_ad, typename T>
InputParameters
FarmsComputeIsotropicElasticityTensorTempl<is_ad, T>::validParams()
{
  InputParameters params = ComputeElasticityTensorBaseTempl<is_ad, T>::validParams();
  params.addClassDescription("Isotropic elasticity tensor with shear modulus updated from initial damage via xi_o and gamma_r.");
  params.addParam<Real>("lambda", 0.0, "Fallback Lamé lambda used when lambda_input is not provided.");
  params.addParam<Real>("shear_modulus", 0.0, "Fallback shear modulus (mu_0) used when shear_modulus_input is not provided.");
  params.addRequiredParam<Real>("xi_o", "Reference strain invariant ratio xi_0 used in gamma_r.");
  params.addParam<std::string>("initial_damage_property", "initial_damage", "Name of the material property holding the initial damage (scalar).");
  params.addParam<MaterialPropertyName>("lambda_input", "", "Optional material property supplying lambda to this object.");
  params.addParam<MaterialPropertyName>("shear_modulus_input", "", "Optional material property supplying shear modulus to this object.");
  return params;
}

template <bool is_ad, typename T>
FarmsComputeIsotropicElasticityTensorTempl<is_ad, T>::FarmsComputeIsotropicElasticityTensorTempl(
    const InputParameters & parameters)
  : ComputeElasticityTensorBaseTempl<is_ad, T>(parameters),
    _lambda_0(this->template getParam<Real>("lambda")),
    _mu_0(this->template getParam<Real>("shear_modulus")),
    _xi_o(this->template getParam<Real>("xi_o")),
    _initial_damage(this->template getMaterialPropertyByName<Real>(
        this->template getParam<std::string>("initial_damage_property"))),
  _lambda_input_prop(nullptr),
  _mu_input_prop(nullptr),
  _use_input_props(false)
{
  // Optional material-property based inputs
  const auto lambda_name = this->template getParam<MaterialPropertyName>("lambda_input");
  const auto mu_name = this->template getParam<MaterialPropertyName>("shear_modulus_input");
  if (!lambda_name.empty() && !mu_name.empty())
  {
    _lambda_input_prop = &this->template getMaterialPropertyByName<Real>(lambda_name);
    _mu_input_prop = &this->template getMaterialPropertyByName<Real>(mu_name);
    _use_input_props = true;
  }
}

template <bool is_ad, typename T>
Real
FarmsComputeIsotropicElasticityTensorTempl<is_ad, T>::computeGammaR(Real lambda0, Real mu0, Real xi0) const
{
  // Same robust formula used elsewhere
  const Real denom = 2.0 * (xi0 * xi0 - 3.0);
  if (std::abs(denom) < 1e-18)
    return 0.0;
  const Real term1 = -xi0 * (-lambda0 * xi0 * xi0 + 6.0 * lambda0 + 2.0 * mu0);
  Real r1 = (lambda0 * xi0 * xi0 + 2.0 * mu0);
  Real r2 = (lambda0 * std::pow(xi0, 4) - 12.0 * lambda0 * xi0 * xi0 + 36.0 * lambda0
             - 6.0 * mu0 * xi0 * xi0 + 24.0 * mu0);
  Real rad = r1 * r2;
  if (rad < 0.0)
    rad = 0.0;
  const Real term2 = std::sqrt(rad);
  const Real gamma_r_minus = (term1 - term2) / denom;
  const Real gamma_r_plus  = (term1 + term2) / denom;
  const Real b_minus = 2.0 * mu0 - gamma_r_minus * xi0;
  return (b_minus > 0.0) ? gamma_r_minus : gamma_r_plus;
}

template <bool is_ad, typename T>
void
FarmsComputeIsotropicElasticityTensorTempl<is_ad, T>::computeQpElasticityTensor()
{
  // Update mu from initial damage at this qp
  const Real d0 = _initial_damage[_qp];
  // Determine base lambda and mu
  Real lambda0 = _lambda_0;
  Real mu0 = _mu_0;
  if (_use_input_props)
  {
    lambda0 = (*_lambda_input_prop)[_qp];
    mu0 = (*_mu_input_prop)[_qp];
  }

  const Real gamma_r = computeGammaR(lambda0, mu0, _xi_o);
  const Real mu_eff = mu0 + _xi_o * d0 * gamma_r;
  const Real lambda = lambda0;

  // Build isotropic tensor from (lambda, mu_eff)
  T C;
  C.fillFromInputVector({lambda, mu_eff}, T::symmetric_isotropic);
  _elasticity_tensor[_qp] = C;

  // Effective stiffness (same approach as isotropic E-nu based estimate using mu)
  _effective_stiffness[_qp] = std::sqrt(std::max(mu_eff, 0.0));
}

// Explicit template instantiation
template class FarmsComputeIsotropicElasticityTensorTempl<false, RankFourTensor>;
template class FarmsComputeIsotropicElasticityTensorTempl<true, RankFourTensor>;
template class FarmsComputeIsotropicElasticityTensorTempl<false, SymmetricRankFourTensor>;
template class FarmsComputeIsotropicElasticityTensorTempl<true, SymmetricRankFourTensor>;
