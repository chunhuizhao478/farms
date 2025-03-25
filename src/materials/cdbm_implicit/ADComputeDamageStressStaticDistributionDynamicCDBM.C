//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADComputeDamageStressStaticDistributionDynamicCDBM.h"

registerMooseObject("farmsApp", ADComputeDamageStressStaticDistributionDynamicCDBM);

InputParameters
ADComputeDamageStressStaticDistributionDynamicCDBM::validParams()
{
  InputParameters params = ADComputeStressBase::validParams();
  params.addClassDescription("Compute stress using elasticity for small strains");
  params.addRequiredParam<Real>("lambda_o","initial lambda value");
  params.addRequiredParam<Real>("shear_modulus_o","initial shear modulus value");
  params.addRequiredParam<Real>("xi_o","xi_o value");
  params.addRequiredParam<Real>("xi_d","xi_d value");
  params.addRequiredParam<Real>("chi","chi value");
  //Compute effective stress
  //------------------------------------------------------------------------------------------------------//
  params.addParam<bool>("compute_effective_stress", false, "Whether or not to compute effective stress");
  params.addParam<Real>("fluid_density", 0, "Fluid density");
  params.addParam<Real>("gravity", 0, "Gravity");
  //------------------------------------------------------------------------------------------------------//
  //Add Cohesion
  //------------------------------------------------------------------------------------------------------//
  params.addParam<bool>("add_cohesion", false, "Whether or not to add cohesion");
  params.addParam<Real>("constant_cohesion", 0, "Constant cohesion");
  params.addParam<Real>("constant_cohesion_cutoff_distance", 0, "Cutoff distance for constant cohesion");
  //------------------------------------------------------------------------------------------------------//
  return params;
}

ADComputeDamageStressStaticDistributionDynamicCDBM::ADComputeDamageStressStaticDistributionDynamicCDBM(const InputParameters & parameters)
  : ADComputeStressBase(parameters),
  _lambda_o(getParam<Real>("lambda_o")),
  _shear_modulus_o(getParam<Real>("shear_modulus_o")),
  _xi_o(getParam<Real>("xi_o")),
  _xi_d(getParam<Real>("xi_d")),
  _chi(getParam<Real>("chi")),
  _initial_damage_val(getADMaterialPropertyByName<Real>("initial_damage")),
  _initial_breakage_val(getADMaterialPropertyByName<Real>("initial_breakage")),
  //Compute effective stress
  //------------------------------------------------------------------------------------------------------//
  _compute_effective_stress(getParam<bool>("compute_effective_stress")),
  _fluid_density(getParam<Real>("fluid_density")),
  _gravity(getParam<Real>("gravity")),
  //------------------------------------------------------------------------------------------------------//
  //Add Cohesion
  //------------------------------------------------------------------------------------------------------//
  _add_cohesion(getParam<bool>("add_cohesion")),
  _constant_cohesion(getParam<Real>("constant_cohesion")),
  _constant_cohesion_cutoff_distance(getParam<Real>("constant_cohesion_cutoff_distance"))
  //------------------------------------------------------------------------------------------------------//
{
  // Check if the parameters are valid
  if (_compute_effective_stress && (_fluid_density == 0 || _gravity == 0))
    mooseError("Fluid density and gravity must be specified when computing effective stress");
  if (_add_cohesion && (_constant_cohesion == 0 || _constant_cohesion_cutoff_distance == 0))
    mooseError("Constant cohesion and its cutoff distance must be specified when adding cohesion");
}

void
ADComputeDamageStressStaticDistributionDynamicCDBM::initialSetup()
{
  // _base_name + "unstabilized_deformation_gradient" is only declared if we're
  // using the Lagrangian kernels.  It's okay to invoke this small strain
  // material if you are using that kernel system and the
  // ComputeLagrangianWrappedStress wrapper
  if (hasBlockMaterialProperty<RankTwoTensor>(_base_name + "strain_increment") &&
      !hasBlockMaterialProperty<RankTwoTensor>(_base_name + "unstabilized_deformation_gradient"))
    mooseError("This linear elastic stress calculation only works for small strains; use "
               "ComputeFiniteStrainElasticStress for simulations using incremental and finite "
               "strains.");
}

void
ADComputeDamageStressStaticDistributionDynamicCDBM::computeQpStress()
{
  
  // Compute gammar
  ADReal gamma_r = computegammar();

  //Compute breakage coefficients
  std::vector<ADReal> avec = computecoefficients(gamma_r);
  ADReal a0 = avec[0];
  ADReal a1 = avec[1];
  ADReal a2 = avec[2];
  ADReal a3 = avec[3];

  // Evaluate shear modulus
  ADReal shear_modulus = _shear_modulus_o + _xi_o * _initial_damage_val[_qp] * gamma_r;
  ADReal gamma_damaged_out = _initial_damage_val[_qp] * gamma_r;

  //
  const ADReal epsilon = 1e-12;
  ADReal I1 = epsilon + _mechanical_strain[_qp](0,0) + _mechanical_strain[_qp](1,1) + _mechanical_strain[_qp](2,2);
  ADReal I2 = epsilon + _mechanical_strain[_qp](0,0) * _mechanical_strain[_qp](0,0) + _mechanical_strain[_qp](1,1) * _mechanical_strain[_qp](1,1) + _mechanical_strain[_qp](2,2) * _mechanical_strain[_qp](2,2) + 2 * _mechanical_strain[_qp](1,2) * _mechanical_strain[_qp](1,2) + 2 * _mechanical_strain[_qp](0,1) * _mechanical_strain[_qp](0,1) + 2 * _mechanical_strain[_qp](0,2) * _mechanical_strain[_qp](0,2);
  ADReal xi = I1 / std::sqrt(I2);

  /* Compute stress */
  ADRankTwoTensor sigma_s = (_lambda_o - gamma_damaged_out / xi) * I1 * ADRankTwoTensor::Identity() + (2 * shear_modulus - gamma_damaged_out * xi) * _mechanical_strain[_qp];
  ADRankTwoTensor sigma_b = (2 * a2 + a1 / xi + 3 * a3 * xi) * I1 * ADRankTwoTensor::Identity() + (2 * a0 + a1 * xi - a3 * std::pow(xi, 3)) * _mechanical_strain[_qp];
  ADRankTwoTensor sigma_total = (1 - _initial_breakage_val[_qp]) * sigma_s + _initial_breakage_val[_qp] * sigma_b;

  // stress
  _stress[_qp] = sigma_total;

  // compute effective stress
  if (_compute_effective_stress){compute_effective_stress();}

  // Add cohesion
  if (_add_cohesion){add_cohesion();}

  // Assign value for elastic strain, which is equal to the mechanical strain
  _elastic_strain[_qp] = _mechanical_strain[_qp];

}

ADReal
ADComputeDamageStressStaticDistributionDynamicCDBM::computegammar()
{
  // Calculate each part of the expression
  ADReal term1 = -_xi_o * (-_lambda_o * std::pow(_xi_o, 2) + 6 * _lambda_o + 2 * _shear_modulus_o);
  ADReal term2_sqrt = std::sqrt((_lambda_o * std::pow(_xi_o, 2) + 2 * _shear_modulus_o) * 
                            (_lambda_o * std::pow(_xi_o, 4) - 12 * _lambda_o * std::pow(_xi_o, 2) + 36 * _lambda_o 
                            - 6 * _shear_modulus_o * std::pow(_xi_o, 2) + 24 * _shear_modulus_o));
  ADReal denominator = 2 * (std::pow(_xi_o, 2) - 3);
  
  // Calculate gamma_r
  ADReal gamma_r = (term1 - term2_sqrt) / denominator;
  
  return gamma_r;
}

std::vector<ADReal>
ADComputeDamageStressStaticDistributionDynamicCDBM::computecoefficients(ADReal gamma_damaged_r)
{

  //compute xi_1
  ADReal _xi_1 = _xi_o + std::sqrt( std::pow(_xi_o , 2) + 2 * _shear_modulus_o / _lambda_o );

  //compute alpha_cr | xi = 0
  ADReal alpha_cr_xi0 = alphacr_root1(0, gamma_damaged_r);

  //compute mu_cr
  ADReal mu_cr = _shear_modulus_o + alpha_cr_xi0 * _xi_o * gamma_damaged_r;

  //a0
  ADReal a0 = _chi * mu_cr;

  //a1
  ADReal numerator_a1 = -2 * _chi * mu_cr * std::pow(_xi_1, 3) + 6 * _chi * mu_cr * _xi_1 * std::pow(_xi_d, 2) - 4 * _chi * mu_cr * std::pow(_xi_d, 3)
                      - 2 * gamma_damaged_r * std::pow(_xi_1, 3) * _xi_d + 2 * gamma_damaged_r * std::pow(_xi_1, 3) * _xi_o
                      + _lambda_o * std::pow(_xi_1, 3) * std::pow(_xi_d, 2) + 2 * _shear_modulus_o * std::pow(_xi_1, 3);
  ADReal denominator_a1 = 2 * std::pow(_xi_1, 3) * _xi_d - 4 * std::pow(_xi_1, 2) * std::pow(_xi_d, 2) + 2 * _xi_1 * std::pow(_xi_d, 3);
  ADReal a1 = numerator_a1 / denominator_a1;

  //a2
  ADReal numerator_a2 = 2 * _chi * mu_cr * std::pow(_xi_1, 3) - 3 * _chi * mu_cr * std::pow(_xi_1, 2) * _xi_d + _chi * mu_cr * std::pow(_xi_d, 3)
                       + 2 * gamma_damaged_r * std::pow(_xi_1, 3) * _xi_d - 2 * gamma_damaged_r * std::pow(_xi_1, 3) * _xi_o
                       - _lambda_o * std::pow(_xi_1, 3) * std::pow(_xi_d, 2) - 2 * _shear_modulus_o * std::pow(_xi_1, 3);
  ADReal denominator_a2 = std::pow(_xi_1, 4) * _xi_d - 2 * std::pow(_xi_1, 3) * std::pow(_xi_d, 2) + std::pow(_xi_1, 2) * std::pow(_xi_d, 3); 
  ADReal a2 = numerator_a2 / denominator_a2; 

  //a3
  ADReal numerator_a3 = -2 * _chi * mu_cr * std::pow(_xi_1, 2) + 4 * _chi * mu_cr * _xi_1 * _xi_d - 2 * _chi * mu_cr * std::pow(_xi_d, 2)
                       - 2 * gamma_damaged_r * std::pow(_xi_1, 2) * _xi_d + 2 * gamma_damaged_r * std::pow(_xi_1, 2) * _xi_o
                       + _lambda_o * std::pow(_xi_1, 2) * std::pow(_xi_d, 2) + 2 * _shear_modulus_o * std::pow(_xi_1, 2);
  ADReal denominator_a3 = 2 * std::pow(_xi_1, 4) * _xi_d - 4 * std::pow(_xi_1, 3) * std::pow(_xi_d, 2) + 2 * std::pow(_xi_1, 2) * std::pow(_xi_d, 3);
  ADReal a3 = numerator_a3 / denominator_a3; 

  //save
  std::vector<ADReal> a_vec {a0,a1,a2,a3};

  return a_vec;

}

// Function for alpha_func_root1
ADReal 
ADComputeDamageStressStaticDistributionDynamicCDBM::alphacr_root1(ADReal xi, ADReal gamma_damaged_r) {
    ADReal term1 = _lambda_o * std::pow(xi, 3) - 6 * _lambda_o * _xi_o + 6 * _shear_modulus_o * xi - 8 * _shear_modulus_o * _xi_o;
    ADReal term2 = std::sqrt(_lambda_o * _lambda_o * std::pow(xi, 6) 
                             - 12 * _lambda_o * _lambda_o * std::pow(xi, 3) * _xi_o 
                             + 36 * _lambda_o * _lambda_o * _xi_o * _xi_o 
                             + 12 * _lambda_o * _shear_modulus_o * std::pow(xi, 4) 
                             - 16 * _lambda_o * _shear_modulus_o * std::pow(xi, 3) * _xi_o 
                             - 72 * _lambda_o * _shear_modulus_o * std::pow(xi, 2) 
                             + 72 * _lambda_o * _shear_modulus_o * xi * _xi_o 
                             + 72 * _lambda_o * _shear_modulus_o 
                             - 12 * _shear_modulus_o * _shear_modulus_o * std::pow(xi, 2) 
                             + 48 * _shear_modulus_o * _shear_modulus_o);
    ADReal denominator = 2 * gamma_damaged_r * (3 * std::pow(xi, 2) - 6 * xi * _xi_o + 4 * _xi_o * _xi_o - 3);
    return (term1 - term2) / denominator;
}

// Function for compute effective stress
void
ADComputeDamageStressStaticDistributionDynamicCDBM::compute_effective_stress()
{
  // Compute effective stress
  ADReal zcoord = _q_point[_qp](2);
  ADReal Pf = _fluid_density * _gravity * std::abs(zcoord);
  _stress[_qp](0,0) += Pf;
  _stress[_qp](1,1) += Pf;
  _stress[_qp](2,2) += Pf;
}

// Function for add cohesion
void
ADComputeDamageStressStaticDistributionDynamicCDBM::add_cohesion()
{
  // Add cohesion
  ADReal zcoord = _q_point[_qp](2);
  if (std::abs(zcoord) < _constant_cohesion_cutoff_distance){
    _stress[_qp](0,0) += _constant_cohesion;
    _stress[_qp](1,1) += _constant_cohesion;
    _stress[_qp](2,2) += _constant_cohesion;
  }
}