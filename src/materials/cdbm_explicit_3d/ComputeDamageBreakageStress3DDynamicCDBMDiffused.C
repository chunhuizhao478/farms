//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ComputeDamageBreakageStress3DDynamicCDBMDiffused.h"
#include "NestedSolve.h"
#include "FEProblem.h"

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
    _gamma_damaged(declareProperty<Real>("damaged_modulus"))
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

  const Real epsilon = 1e-12;
  Real I1 = epsilon + eps_e(0,0) + eps_e(1,1) + eps_e(2,2);
  Real I2 = epsilon + eps_e(0,0) * eps_e(0,0) + eps_e(1,1) * eps_e(1,1) + eps_e(2,2) * eps_e(2,2) + 2 * eps_e(0,1) * eps_e(0,1) + 2 * eps_e(0,2) * eps_e(0,2) + 2 * eps_e(1,2) * eps_e(1,2);
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
  computeQpTangentModulus(tangent,I1,I2,xi,eps_e);
  _Jacobian_mult[_qp] = tangent;

  // Compute deviatoric strain rate
  computeDeviatroicStrainRateTensor();

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
                                                      RankTwoTensor Ee)
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

  const Real adjusted_I2 = (I2 <= 1e-12) ? 1e-12 : I2;
  const RankTwoTensor identity = RankTwoTensor::Identity();

  // Precompute dxidE tensor
  RankTwoTensor dxidE_tensor;
  for (unsigned int k = 0; k < 3; ++k)
    for (unsigned int l = 0; l < 3; ++l)
      dxidE_tensor(k, l) = (identity(k, l) * adjusted_I2 - I1 * Ee(k, l)) / std::pow(adjusted_I2, 1.5);

  const RankTwoTensor dxim1dE_tensor = dxidE_tensor * (-1.0 / (xi * xi));

  // Compute terms for dSedE
  const Real lambda_term = _lambda[_qp] - _gamma_damaged[_qp] / xi;
  const Real shear_term = 2.0 * _shear_modulus[_qp] - _gamma_damaged[_qp] * xi;

  RankFourTensor term_se1 = identity.outerProduct(-_gamma_damaged[_qp] * I1 * dxim1dE_tensor);
  RankFourTensor term_se2 = identity.outerProduct(identity) * lambda_term;
  RankFourTensor term_se3 = Ee.outerProduct(-_gamma_damaged[_qp] * dxidE_tensor);
  RankFourTensor term_se4 = RankFourTensor(RankFourTensor::initIdentityFour) * shear_term;

  RankFourTensor dSedE = term_se1 + term_se2 + term_se3 + term_se4;

  // Compute terms for dSbdE
  const Real coeff2_b = 2.0 * a2 + a1 / xi + 3.0 * a3 * xi;
  const Real coeff4_b = 2.0 * a0 + a1 * xi - a3 * xi * xi * xi;

  RankFourTensor term_b1 = identity.outerProduct((a1 * dxim1dE_tensor + 3 * a3 * dxidE_tensor) * I1);
  RankFourTensor term_b2 = identity.outerProduct(identity) * coeff2_b;
  RankFourTensor term_b3 = Ee.outerProduct(a1 * dxidE_tensor - a3 * 3 * xi * xi * dxidE_tensor);
  RankFourTensor term_b4 = RankFourTensor(RankFourTensor::initIdentityFour) * coeff4_b;

  RankFourTensor dSbdE = term_b1 + term_b2 + term_b3 + term_b4;

  // Combine and assign tangent
  tangent = (1.0 - _B_damagedvar_aux[_qp]) * dSedE + _B_damagedvar_aux[_qp] * dSbdE;  

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