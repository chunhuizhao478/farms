//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ComputeDamageBreakageStress3DDynamicCDBMDiffused.h"

registerMooseObject("farmsApp", ComputeDamageBreakageStress3DDynamicCDBMDiffused);

InputParameters
ComputeDamageBreakageStress3DDynamicCDBMDiffused::validParams()
{ 
  //Note: lambda_o, shear_modulus_o is defined in "ComputeGeneralDamageBreakageStressBase"
  //to initialize _lambda, _shear_modulus material properties
  InputParameters params = ComputeDamageBreakageStressBase3D::validParams();
  params.addClassDescription("Compute stress using elasticity for small strains");
  
  //constant parameters
  params.addRequiredParam<Real>(        "lambda_o", "initial lambda constant value");
  params.addRequiredParam<Real>( "shear_modulus_o", "initial shear modulus value");
  params.addRequiredParam<Real>(            "xi_0", "strain invariants ratio: onset of damage evolution");
  params.addRequiredParam<Real>(            "xi_d", "strain invariants ratio: onset of breakage healing");
  params.addRequiredParam<Real>(             "chi", "ratio of solid energy and granular energy");
  params.addRequiredParam<Real>(             "C_g", "material parameter: compliance or fluidity of the fine grain granular material");
  params.addRequiredParam<Real>(              "m1", "coefficient of std::power law indexes");
  params.addRequiredParam<Real>(              "m2", "coefficient of std::power law indexes");

  //input coupled variables from sub app
  params.addRequiredCoupledVar("alpha_damagedvar_aux", "second_elastic_strain_invariant");
  params.addRequiredCoupledVar("B_damagedvar_aux", "strain_invariant_ratio");

  //porous flow coupling
  params.addParam<bool>("porous_flow_coupling", false, "Flag to indicate if this material is coupled to a porous flow module");
  params.addParam<Real>("coeff_b", -1.0, "the factor to govern the exponential permeability enhancement");
  params.addParam<Real>("intrinsic_permeability", -1.0, "the intrisic permeability for the material, the intrisic permeability is used to take diagonal of the permeability tensor");

  return params;
}

ComputeDamageBreakageStress3DDynamicCDBMDiffused::ComputeDamageBreakageStress3DDynamicCDBMDiffused(const InputParameters & parameters)
  : ComputeStressBase(parameters),
    _lambda_o(getParam<Real>("lambda_o")),
    _shear_modulus_o(getParam<Real>("shear_modulus_o")),
    _xi_0(getParam<Real>("xi_0")),
    _xi_d(getParam<Real>("xi_d")),
    _chi(getParam<Real>("chi")),
    _C_g(getParam<Real>("C_g")),
    _m1(getParam<Real>("m1")),
    _m2(getParam<Real>("m2")),
    _alpha_damagedvar_aux(coupledValue("alpha_damagedvar_aux")),
    _B_damagedvar_aux(coupledValue("B_damagedvar_aux")),
    _deviatroic_strain_rate(declareProperty<Real>("deviatroic_strain_rate")),
    _sigma_d(declareProperty<RankTwoTensor>("deviatoric_stress")),
    _sigma_d_old(getMaterialPropertyOldByName<RankTwoTensor>("deviatoric_stress")),
    _eps_p(declareProperty<RankTwoTensor>("plastic_strain_tensor")),
    _eps_p_old(getMaterialPropertyOldByName<RankTwoTensor>("plastic_strain_tensor")),
    _eps_total(declareProperty<RankTwoTensor>("total_strain_tensor")),
    _eps_total_old(getMaterialPropertyOldByName<RankTwoTensor>("total_strain_tensor")),
    _eps_e(declareProperty<RankTwoTensor>("elastic_strain_tensor")),
    _I1(declareProperty<Real>("first_elastic_strain_invariant")),
    _I2(declareProperty<Real>("second_elastic_strain_invariant")),
    _xi(declareProperty<Real>("strain_invariant_ratio")),
    _lambda(declareProperty<Real>("lambda_const")),
    _shear_modulus(declareProperty<Real>("shear_modulus")),
    _gamma_damaged(declareProperty<Real>("damaged_modulus")),
    _stress_perturbation(getMaterialPropertyOldByName<Real>("shear_stress_perturbation")),
    //## Porous flow coupling ##
    _porous_flow_coupling(getParam<bool>("porous_flow_coupling")),
    _solid_bulk_compliance_damaged(declareProperty<Real>("solid_bulk_compliance_damaged")),
    _coeff_b(getParam<Real>("coeff_b")),
    _intrinsic_permeability(getParam<Real>("intrinsic_permeability")),
    _effective_perm(declareProperty<RealTensorValue>("effective_perm")),
    _crack_rotation(declareProperty<RankTwoTensor>("crack_rotation")),
    _elasticity_tensor_name(_base_name + "elasticity_tensor"),
    _elasticity_tensor(getMaterialPropertyByName<RankFourTensor>(_elasticity_tensor_name))
{
}

void
ComputeDamageBreakageStress3DDynamicCDBMDiffused::initialSetup()
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
ComputeDamageBreakageStress3DDynamicCDBMDiffused::initQpStatefulProperties()
{
  _elastic_strain[_qp].zero();
  _stress[_qp].zero();
}

void
ComputeDamageBreakageStress3DDynamicCDBMDiffused::computeQpStress()
{ 
  
  /*
  compute gammar, breakage coefficients
  */
  Real gamma_damaged_r = computegammar();
  std::vector<Real> avec = computecoefficients(gamma_damaged_r);
  Real a0 = avec[0];
  Real a1 = avec[1];
  Real a2 = avec[2];
  Real a3 = avec[3];

  //lambda, shear_modulus, gamma_damaged are updated
  Real lambda_out = _lambda_o;
  Real shear_modulus_out = _shear_modulus_o + _alpha_damagedvar_aux[_qp] * _xi_0 * gamma_damaged_r;
  Real gamma_damaged_out = _alpha_damagedvar_aux[_qp] * gamma_damaged_r;

  _lambda[_qp] = lambda_out;
  _shear_modulus[_qp] = shear_modulus_out;
  _gamma_damaged[_qp] = gamma_damaged_out;

  /* compute strain */
  RankTwoTensor eps_p = _eps_p_old[_qp] + _dt * _C_g * std::pow(_B_damagedvar_aux[_qp],_m1) * _sigma_d_old[_qp];
  RankTwoTensor eps_e = _mechanical_strain[_qp] - eps_p;

  Real I1 = eps_e(0,0) + eps_e(1,1) + eps_e(2,2);
  Real I2 = eps_e(0,0) * eps_e(0,0) + eps_e(1,1) * eps_e(1,1) + eps_e(2,2) * eps_e(2,2) + 2 * eps_e(0,1) * eps_e(0,1) + 2 * eps_e(0,2) * eps_e(0,2) + 2 * eps_e(1,2) * eps_e(1,2);
  Real xi = I1/std::sqrt(I2); //catch the nan error in the initial solve

  //Represent sigma (solid(s) + granular(b))
  RankTwoTensor sigma_s;
  RankTwoTensor sigma_b;
  RankTwoTensor sigma_total;
  RankTwoTensor sigma_d;
  const auto I = RankTwoTensor::Identity();

  /* Compute stress */
  sigma_s = (lambda_out - gamma_damaged_out / xi) * I1 * RankTwoTensor::Identity() + (2 * shear_modulus_out - gamma_damaged_out * xi) * eps_e;
  sigma_b = (2 * a2 + a1 / xi + 3 * a3 * xi) * I1 * RankTwoTensor::Identity() + (2 * a0 + a1 * xi - a3 * std::pow(xi, 3)) * eps_e;
  sigma_total = (1 - _B_damagedvar_aux[_qp]) * sigma_s + _B_damagedvar_aux[_qp] * sigma_b;

  /* add pore pressure */
  if (_stress_perturbation[_qp] != 0){
    sigma_total(0,0) -= _stress_perturbation[_qp];
    sigma_total(1,1) -= _stress_perturbation[_qp];
    sigma_total(2,2) -= _stress_perturbation[_qp];
  }

  sigma_d = sigma_total - 0.3333 * (sigma_total(0,0) + sigma_total(1,1) + sigma_total(2,2)) * I;

  _eps_total[_qp] = eps_p + eps_e;
  _eps_p[_qp] = eps_p;
  _eps_e[_qp] = eps_e;
  _I1[_qp] = I1;
  _I2[_qp] = I2;
  _xi[_qp] = xi;
  _sigma_d[_qp] = sigma_d;

  // Rotate the stress state to the current configuration
  _stress[_qp] = sigma_total;

  // Assign value for elastic strain, which is equal to the mechanical strain
  _elastic_strain[_qp] = eps_e;

  // Compute tangent
  RankFourTensor tangent;
  computeQpTangentModulus(tangent,I1,I2,xi,eps_e,a0,a1,a2,a3,gamma_damaged_r);
  _Jacobian_mult[_qp] = tangent;

  // Compute deviatoric strain rate
  computeDeviatroicStrainRateTensor();

  /* Compute Principal Strains and Rotation Matrix */
  RealVectorValue strain_in_crack_dir; //principal strains
  computeCrackStrainAndOrientation(strain_in_crack_dir);

  // Compute bulk compliance
  updateSolidBulkCompliance();

  // Compute effective permeability
  updatePermeabilityForCracking();

}

Real 
ComputeDamageBreakageStress3DDynamicCDBMDiffused::computegammar()
{
  // Calculate each part of the expression
  Real term1 = -_xi_0 * (-_lambda_o * pow(_xi_0, 2) + 6 * _lambda_o + 2 * _shear_modulus_o);
  Real term2_sqrt = sqrt((_lambda_o * pow(_xi_0, 2) + 2 * _shear_modulus_o) * 
                            (_lambda_o * pow(_xi_0, 4) - 12 * _lambda_o * pow(_xi_0, 2) + 36 * _lambda_o 
                            - 6 * _shear_modulus_o * pow(_xi_0, 2) + 24 * _shear_modulus_o));
  Real denominator = 2 * (pow(_xi_0, 2) - 3);
  
  // Calculate gamma_r
  Real gamma_r = (term1 - term2_sqrt) / denominator;
  
  //save
  return gamma_r;
}

std::vector<Real>
ComputeDamageBreakageStress3DDynamicCDBMDiffused::computecoefficients(Real gamma_damaged_r)
{

  //compute xi_1
  Real _xi_1 = _xi_0 + sqrt( pow(_xi_0 , 2) + 2 * _shear_modulus_o / _lambda_o );

  //compute alpha_cr | xi = 0
  Real alpha_cr_xi0 = alphacr_root1(0, gamma_damaged_r);

  //compute mu_cr
  Real mu_cr = _shear_modulus_o + alpha_cr_xi0 * _xi_0 * gamma_damaged_r;

  //a0
  Real a0 = _chi * mu_cr;

  //a1
  Real numerator_a1 = -2 * _chi * mu_cr * pow(_xi_1, 3) + 6 * _chi * mu_cr * _xi_1 * pow(_xi_d, 2) - 4 * _chi * mu_cr * pow(_xi_d, 3)
                      - 2 * gamma_damaged_r * pow(_xi_1, 3) * _xi_d + 2 * gamma_damaged_r * pow(_xi_1, 3) * _xi_0
                      + _lambda_o * pow(_xi_1, 3) * pow(_xi_d, 2) + 2 * _shear_modulus_o * pow(_xi_1, 3);
  Real denominator_a1 = 2 * pow(_xi_1, 3) * _xi_d - 4 * pow(_xi_1, 2) * pow(_xi_d, 2) + 2 * _xi_1 * pow(_xi_d, 3);
  Real a1 = numerator_a1 / denominator_a1;

  //a2
  Real numerator_a2 = 2 * _chi * mu_cr * pow(_xi_1, 3) - 3 * _chi * mu_cr * pow(_xi_1, 2) * _xi_d + _chi * mu_cr * pow(_xi_d, 3)
                       + 2 * gamma_damaged_r * pow(_xi_1, 3) * _xi_d - 2 * gamma_damaged_r * pow(_xi_1, 3) * _xi_0
                       - _lambda_o * pow(_xi_1, 3) * pow(_xi_d, 2) - 2 * _shear_modulus_o * pow(_xi_1, 3);
  Real denominator_a2 = pow(_xi_1, 4) * _xi_d - 2 * pow(_xi_1, 3) * pow(_xi_d, 2) + pow(_xi_1, 2) * pow(_xi_d, 3); 
  Real a2 = numerator_a2 / denominator_a2; 

  //a3
  Real numerator_a3 = -2 * _chi * mu_cr * pow(_xi_1, 2) + 4 * _chi * mu_cr * _xi_1 * _xi_d - 2 * _chi * mu_cr * pow(_xi_d, 2)
                       - 2 * gamma_damaged_r * pow(_xi_1, 2) * _xi_d + 2 * gamma_damaged_r * pow(_xi_1, 2) * _xi_0
                       + _lambda_o * pow(_xi_1, 2) * pow(_xi_d, 2) + 2 * _shear_modulus_o * pow(_xi_1, 2);
  Real denominator_a3 = 2 * pow(_xi_1, 4) * _xi_d - 4 * pow(_xi_1, 3) * pow(_xi_d, 2) + 2 * pow(_xi_1, 2) * pow(_xi_d, 3);
  Real a3 = numerator_a3 / denominator_a3; 

  //save
  std::vector<Real> a_vec {a0,a1,a2,a3};

  return a_vec;

}

// Function for alpha_func_root1
Real 
ComputeDamageBreakageStress3DDynamicCDBMDiffused::alphacr_root1(Real xi, Real gamma_damaged_r) {
    Real term1 = _lambda_o * pow(xi, 3) - 6 * _lambda_o * _xi_0 + 6 * _shear_modulus_o * xi - 8 * _shear_modulus_o * _xi_0;
    Real term2 = std::sqrt(_lambda_o * _lambda_o * pow(xi, 6) 
                             - 12 * _lambda_o * _lambda_o * pow(xi, 3) * _xi_0 
                             + 36 * _lambda_o * _lambda_o * _xi_0 * _xi_0 
                             + 12 * _lambda_o * _shear_modulus_o * pow(xi, 4) 
                             - 16 * _lambda_o * _shear_modulus_o * pow(xi, 3) * _xi_0 
                             - 72 * _lambda_o * _shear_modulus_o * pow(xi, 2) 
                             + 72 * _lambda_o * _shear_modulus_o * xi * _xi_0 
                             + 72 * _lambda_o * _shear_modulus_o 
                             - 12 * _shear_modulus_o * _shear_modulus_o * pow(xi, 2) 
                             + 48 * _shear_modulus_o * _shear_modulus_o);
    Real denominator = 2 * gamma_damaged_r * (3 * pow(xi, 2) - 6 * xi * _xi_0 + 4 * _xi_0 * _xi_0 - 3);
    return (term1 - term2) / denominator;
}

// Function for alpha_func_root2
Real 
ComputeDamageBreakageStress3DDynamicCDBMDiffused::alphacr_root2(Real xi, Real gamma_damaged_r) {
    return 2 * _shear_modulus_o / (gamma_damaged_r * (xi - 2 * _xi_0));
}

void
ComputeDamageBreakageStress3DDynamicCDBMDiffused::computeQpTangentModulus(RankFourTensor & tangent, 
                                                      Real I1, 
                                                      Real I2, 
                                                      Real xi, 
                                                      RankTwoTensor Ee,
                                                      Real a0,
                                                      Real a1,
                                                      Real a2,
                                                      Real a3,
                                                      Real gamma_damaged_r)
{

  // Use consistent values - same as in stress computation

  // Use the SAME values as in stress computation
  Real lambda_out = _lambda_o;
  Real shear_modulus_out = _shear_modulus_o + _alpha_damagedvar_aux[_qp] * _xi_0 * gamma_damaged_r;
  Real gamma_damaged_out = _alpha_damagedvar_aux[_qp] * gamma_damaged_r;

  // Safety check for small I2
  const Real adjusted_I2 = std::max(I2, 1e-12);
  const Real sqrt_I2 = std::sqrt(adjusted_I2);
  const RankTwoTensor identity = RankTwoTensor::Identity();

  // Check for limiting case: alpha = 0, B = 0 (should return elasticity tensor)
  // if (std::abs(_alpha_damagedvar_aux[_qp]) < 1e-12 && std::abs(_B_damagedvar_aux[_qp]) < 1e-12) {
  //   // Standard elasticity tensor: C_ijkl = λ δ_ij δ_kl + μ (δ_ik δ_jl + δ_il δ_jk)
  //   tangent.zero();
  //   for (unsigned int i = 0; i < 3; ++i) {
  //     for (unsigned int j = 0; j < 3; ++j) {
  //       for (unsigned int k = 0; k < 3; ++k) {
  //         for (unsigned int l = 0; l < 3; ++l) {
  //           tangent(i, j, k, l) = _lambda_o * identity(i, j) * identity(k, l) + 
  //                                 _shear_modulus_o * (identity(i, k) * identity(j, l) + identity(i, l) * identity(j, k));
  //         }
  //       }
  //     }
  //   }
  //   return;
  // }

  // Corrected derivative: ∂ξ/∂E_kl = ∂(I1/√I2)/∂E_kl
  // = (∂I1/∂E_kl * √I2 - I1 * ∂I2/∂E_kl / (2√I2)) / I2
  // where ∂I1/∂E_kl = δ_kl and ∂I2/∂E_kl = 2*E_kl
  RankTwoTensor dxidE_tensor;
  for (unsigned int k = 0; k < 3; ++k) {
    for (unsigned int l = 0; l < 3; ++l) {
      dxidE_tensor(k, l) = identity(k, l) / sqrt_I2 - I1 * Ee(k, l) / std::pow(adjusted_I2, 1.5);
    }
  }

  // ∂(1/ξ)/∂E = -1/ξ² * ∂ξ/∂E
  const RankTwoTensor dxim1dE_tensor = dxidE_tensor * (-1.0 / (xi * xi));

  // Compute solid phase tangent (dSs/dE)
  const Real lambda_term = lambda_out - gamma_damaged_out / xi;
  const Real shear_term = 2.0 * shear_modulus_out - gamma_damaged_out * xi;

  RankFourTensor dSsdE;
  dSsdE.zero();
  
  // CORRECTED: Complete implementation of solid phase tangent
  // ∂S^s_ij/∂E_kl = (-γ ∂ξ^(-1)/∂E_kl)I_1 δ_ij + (λ - γ/ξ) ∂I_1/∂E_kl δ_ij + (-γ ∂ξ/∂E_kl)E_ij + (2μ - γξ) ∂E_ij/∂E_kl
  for (unsigned int i = 0; i < 3; ++i) {
    for (unsigned int j = 0; j < 3; ++j) {
      for (unsigned int k = 0; k < 3; ++k) {
        for (unsigned int l = 0; l < 3; ++l) {
          // Term 1a: (λ - γ/ξ) * ∂I1/∂E_kl * δ_ij = (λ - γ/ξ) * δ_kl * δ_ij
          dSsdE(i, j, k, l) += lambda_term * identity(i, j) * identity(k, l);
          
          // Term 1b: (-γ ∂ξ^(-1)/∂E_kl) * I1 * δ_ij - PREVIOUSLY MISSING
          dSsdE(i, j, k, l) -= gamma_damaged_out * dxim1dE_tensor(k, l) * I1 * identity(i, j);
          
          // Term 2a: (2μ - γξ) * ∂E_ij/∂E_kl
          Real I4_ijkl = 0.5 * (identity(i, k) * identity(j, l) + identity(i, l) * identity(j, k));
          dSsdE(i, j, k, l) += shear_term * I4_ijkl;
          
          // Term 2b: (-γ ∂ξ/∂E_kl) * E_ij
          dSsdE(i, j, k, l) -= gamma_damaged_out * dxidE_tensor(k, l) * Ee(i, j);
        }
      }
    }
  }

  // Compute granular phase tangent (dSb/dE)
  const Real coeff2_b = 2.0 * a2 + a1 / xi + 3.0 * a3 * xi;
  const Real coeff4_b = 2.0 * a0 + a1 * xi - a3 * xi * xi * xi;

  RankFourTensor dSbdE;
  dSbdE.zero();
  
 // CORRECTED: Complete implementation of granular phase tangent
  // ∂S^b_ij/∂E_kl = (a_1 ∂ξ^(-1)/∂E_kl + 3a_3 ∂ξ/∂E_kl)I_1 δ_ij + (2a_2 + a_1/ξ + 3a_3ξ) ∂I_1/∂E_kl δ_ij
  //                + (a_1 ∂ξ/∂E_kl - a_3 ∂ξ^3/∂E_kl)E_ij + (2a_0 + a_1ξ - a_3ξ^3) ∂E_ij/∂E_kl
  for (unsigned int i = 0; i < 3; ++i) {
    for (unsigned int j = 0; j < 3; ++j) {
      for (unsigned int k = 0; k < 3; ++k) {
        for (unsigned int l = 0; l < 3; ++l) {
          // Term 1a: (2a_2 + a_1/ξ + 3a_3ξ) * ∂I1/∂E_kl * δ_ij = coeff2_b * δ_kl * δ_ij
          dSbdE(i, j, k, l) += coeff2_b * identity(i, j) * identity(k, l);
          
          // Term 1b: a_1 * ∂ξ^(-1)/∂E_kl * I1 * δ_ij - PREVIOUSLY MISSING
          dSbdE(i, j, k, l) += a1 * dxim1dE_tensor(k, l) * I1 * identity(i, j);
          
          // Term 1c: 3a_3 * ∂ξ/∂E_kl * I1 * δ_ij
          dSbdE(i, j, k, l) += 3.0 * a3 * dxidE_tensor(k, l) * I1 * identity(i, j);
          
          // Term 2a: (2a_0 + a_1ξ - a_3ξ^3) * ∂E_ij/∂E_kl
          Real I4_ijkl = 0.5 * (identity(i, k) * identity(j, l) + identity(i, l) * identity(j, k));
          dSbdE(i, j, k, l) += coeff4_b * I4_ijkl;
          
          // Term 2b: a_1 * ∂ξ/∂E_kl * E_ij
          dSbdE(i, j, k, l) += a1 * dxidE_tensor(k, l) * Ee(i, j);
          
          // Term 2c: -a_3 * ∂ξ^3/∂E_kl * E_ij
          // ∂ξ^3/∂E_kl = 3ξ^2 * ∂ξ/∂E_kl
          dSbdE(i, j, k, l) -= a3 * 3.0 * xi * xi * dxidE_tensor(k, l) * Ee(i, j);
        }
      }
    }
  }

  // Combine: tangent = (1-B)*dSs/dE + B*dSb/dE
  tangent = dSsdE * (1.0 - _B_damagedvar_aux[_qp]) + dSbdE * _B_damagedvar_aux[_qp]; 

}

void
ComputeDamageBreakageStress3DDynamicCDBMDiffused::computeDeviatroicStrainRateTensor()
{
  //Compute strain rate E_dot = F^T * D * F
  RankTwoTensor E_dot = (_eps_total[_qp] - _eps_total_old[_qp]) / _dt;
  //Compute deviatoric strain rate tensor E_dev_dot 
  RankTwoTensor E_dev_dot = E_dot - (1.0/3.0) * E_dot.trace() * RankTwoTensor::Identity();
  //Compute J2_dot = 1/2 * E_dev_dot(i,j) * E_dev_dot(i,j)
  Real J2_dot = 0.0;
  for (unsigned int i = 0; i < 3; ++i){
    for (unsigned int j = 0; j < 3; ++j){
      J2_dot += 0.5 * E_dev_dot(i,j) * E_dev_dot(i,j);
    }
  }
  //Compute equivalent strain rate
  _deviatroic_strain_rate[_qp] = std::sqrt(2.0/3.0 * J2_dot);
}

void
ComputeDamageBreakageStress3DDynamicCDBMDiffused::computeCrackStrainAndOrientation(
    RealVectorValue & strain_in_crack_dir)
{
  // The rotation tensor is ordered such that directions for pre-existing cracks appear first
  // in the list of columns.  For example, if there is one existing crack, its direction is in the
  // first column in the rotation tensor.
  
  // If porous flow coupling is not enabled, return
  if (!_porous_flow_coupling)
  return;

  std::vector<Real> eigval(3, 0.0);
  RankTwoTensor eigvec;

  _elastic_strain[_qp].symmetricEigenvaluesEigenvectors(eigval, eigvec);

  // If the elastic strain is beyond the cracking strain, save the eigen vectors as
  // the rotation tensor. Reverse their order so that the third principal strain
  // (most tensile) will correspond to the first crack.
  _crack_rotation[_qp].fillColumn(0, eigvec.column(2));
  _crack_rotation[_qp].fillColumn(1, eigvec.column(1));
  _crack_rotation[_qp].fillColumn(2, eigvec.column(0));

  strain_in_crack_dir(0) = eigval[2];
  strain_in_crack_dir(1) = eigval[1];
  strain_in_crack_dir(2) = eigval[0];
}

void
ComputeDamageBreakageStress3DDynamicCDBMDiffused::updateSolidBulkCompliance()
{

  // If porous flow coupling is not enabled, return
  if (!_porous_flow_coupling)
  return;

  /*
  compute gammar, breakage coefficients
  */
  Real gamma_damaged_r = computegammar();
  std::vector<Real> avec = computecoefficients(gamma_damaged_r);
  Real a0 = avec[0];
  Real a1 = avec[1];
  Real a2 = avec[2];
  Real a3 = avec[3];

  // bulk modulus 
  // solid phase bulk modulus: K_s = lambda_eff + 2/3 * shear_modulus_eff = (lambda - gamma_damaged / xi) + 2/3 * (shear_modulus - gamma_damaged * xi / 2)
  // granular phase bulk modulus: K_b = (2 * a2 + a1 / xi + 3 * a3 * xi) + 2/3 * (2 * a0 + a1 * xi - a3 * xi^3)
  // effective bulk modulus: K_eff = (1 - B_damagedvar_aux) * K_s + B_damagedvar_aux * K_b

  Real K_s = (_lambda[_qp] - _gamma_damaged[_qp] / _xi[_qp]) + 2.0 / 3.0 * (_shear_modulus[_qp] - _gamma_damaged[_qp] * _xi[_qp] / 2.0);
  Real K_b = (2.0 * a2 + a1 / _xi[_qp] + 3.0 * a3 * _xi[_qp]) + 2.0 / 3.0 * (2.0 * a0 + a1 * _xi[_qp] - a3 * std::pow(_xi[_qp], 3));
  Real K_eff = (1.0 - _B_damagedvar_aux[_qp]) * K_s + _B_damagedvar_aux[_qp] * K_b;

  //compute compliance
  _solid_bulk_compliance_damaged[_qp] = 1.0 / K_eff;

}

void
ComputeDamageBreakageStress3DDynamicCDBMDiffused::updatePermeabilityForCracking()
{

  // If porous flow coupling is not enabled, return
  if (!_porous_flow_coupling)
    return;

  // Get transformation matrix
  const RankTwoTensor & R = _crack_rotation[_qp];

  // Initialize effective permeability new
  RankTwoTensor effective_perm_new;

  //Compute the intrinsic permeability
  RankTwoTensor perm_intrinsic = _intrinsic_permeability * RankTwoTensor::Identity();

  // Initialize effective permeability new
  effective_perm_new = perm_intrinsic * std::exp( _alpha_damagedvar_aux[_qp] * _coeff_b );

  // Rotate back to global frame
  effective_perm_new.rotate(R);

  // Update effective perm
  _effective_perm[_qp] = effective_perm_new;

}

void
ComputeDamageBreakageStress3DDynamicCDBMDiffused::HeavisideFunction(Real & x)
{
  if (x > 0.0)
    x = 1.0;
  else
    x = 0.0;
}