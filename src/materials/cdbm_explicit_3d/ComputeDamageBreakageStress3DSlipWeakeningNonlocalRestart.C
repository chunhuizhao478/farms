//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ComputeDamageBreakageStress3DSlipWeakeningNonlocalRestart.h"
#include "NestedSolve.h"
#include "FEProblem.h"

registerMooseObject("farmsApp", ComputeDamageBreakageStress3DSlipWeakeningNonlocalRestart);

InputParameters
ComputeDamageBreakageStress3DSlipWeakeningNonlocalRestart::validParams()
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
  params.addRequiredParam<Real>(          "xi_min", "strain invariants ratio: minimum allowable value");
  params.addRequiredParam<Real>(          "xi_max", "strain invariants ratio: maximum allowable value");
  params.addRequiredParam<Real>(             "chi", "ratio of solid energy and granular energy");
  params.addRequiredParam<Real>(             "C_g", "material parameter: compliance or fluidity of the fine grain granular material");
  params.addRequiredParam<Real>(              "m1", "coefficient of std::power law indexes");
  params.addRequiredParam<Real>(              "m2", "coefficient of std::power law indexes");
  params.addRequiredParam<Real>(     "Cd_constant", "coefficient gives positive damage evolution");
  params.addRequiredParam<Real>(             "C_1", "coefficient of healing for damage evolution");
  params.addRequiredParam<Real>(             "C_2", "coefficient of healing for damage evolution");
  params.addRequiredParam<Real>(      "beta_width", "coefficient gives width of transitional region");
  params.addRequiredParam<Real>( "CdCb_multiplier", "multiplier between Cd and Cb");
  params.addRequiredParam<Real>(    "CBH_constant", "constant CBH value");

  //strain rate dependent Cd parameters
  params.addParam<bool>("use_strain_rate_dependent_Cd", false,
                        "Use strain rate dependent Cd (default: false)");
  params.addParam<Real>( "m_exponent", 0.8, "strain rate dependent parameters");
  params.addParam<Real>( "strain_rate_hat", 1e-4, "strain rate dependent parameters");
  params.addParam<Real>( "cd_hat", 1.0, "strain rate dependent parameters");
  params.addParam<bool>("zero_Cd_below_threshold", false,
                        "If true, set Cd = 0 when deviatoric strain rate < strain_rate_hat; otherwise use cd_hat (default: false).");

  //use nonlocal equivalent strain
  params.addParam<bool>("use_nonlocal_eqstrain", false,
                        "Use nonlocal equivalent strain (default: false)");
  params.addParam<std::vector<unsigned int>>("nonlocal_eqstrain_blocks", {},
                        "REQUIRED when use_nonlocal_eqstrain=true. Subdomain/Block IDs where nonlocal equivalent strain is enabled (e.g., 100 200)");

  //use nonlocal strain rate for Cd calculation
  params.addParam<bool>("use_nonlocal_strain_rate", false,
                        "Use nonlocal averaged strain rate for Cd calculation (default: false)");

  //static solve flag
  params.addParam<bool>("static_solve_flag", true,
                        "Flag to determine which part of setupInitial() to use (default: true)");

  return params;
}

ComputeDamageBreakageStress3DSlipWeakeningNonlocalRestart::ComputeDamageBreakageStress3DSlipWeakeningNonlocalRestart(const InputParameters & parameters)
  : ComputeDamageBreakageStressBase3D(parameters),
    _xi_0(getParam<Real>("xi_0")),
    _xi_d(getParam<Real>("xi_d")),
    _xi_min(getParam<Real>("xi_min")),
    _xi_max(getParam<Real>("xi_max")),
    _chi(getParam<Real>("chi")),
    _C_g(getParam<Real>("C_g")),
    _m1(getParam<Real>("m1")),
    _m2(getParam<Real>("m2")),
    _alpha_damagedvar_old(getMaterialPropertyOldByName<Real>("alpha_damagedvar")),
    _B_old(getMaterialPropertyOldByName<Real>("B")),
    _xi_old(getMaterialPropertyOldByName<Real>("xi")),
    _I1_old(getMaterialPropertyOldByName<Real>("I1")),
    _I2_old(getMaterialPropertyOldByName<Real>("I2")),
    _lambda_old(getMaterialPropertyOldByName<Real>("lambda")),
    _shear_modulus_old(getMaterialPropertyOldByName<Real>("shear_modulus")),
    _gamma_damaged_old(getMaterialPropertyOldByName<Real>("gamma_damaged")),
    _eps_total_old(getMaterialPropertyOldByName<RankTwoTensor>("eps_total")),
    _mechanical_strain_old(getMaterialPropertyOldByName<RankTwoTensor>("mechanical_strain")),
    _eps_p_old(getMaterialPropertyOldByName<RankTwoTensor>("eps_p")),
    _eps_e_old(getMaterialPropertyOldByName<RankTwoTensor>("eps_e")),
    _sigma_d_old(getMaterialPropertyOldByName<RankTwoTensor>("sigma_d")),
    _sts_total_old(getMaterialPropertyOldByName<RankTwoTensor>("sts_total")),
    _static_initial_stress_tensor(getMaterialProperty<RankTwoTensor>("static_initial_stress_tensor")),
    _static_initial_strain_tensor(getMaterialProperty<RankTwoTensor>("static_initial_strain_tensor")),
    _sts_initial_tensor_old(getMaterialPropertyOldByName<RankTwoTensor>("sts_initial_tensor")),
    _initial_damage(getMaterialPropertyByName<Real>("initial_damage")),
    _initial_breakage(getMaterialPropertyByName<Real>("initial_breakage")),
    _damage_perturbation(getMaterialPropertyByName<Real>("damage_perturbation")),
    _Cd_constant(getParam<Real>("Cd_constant")),
    _C1(getParam<Real>("C_1")),
    _C2(getParam<Real>("C_2")),
    _beta_width(getParam<Real>("beta_width")),
    _CdCb_multiplier(getParam<Real>("CdCb_multiplier")),
    _CBH_constant(getParam<Real>("CBH_constant")),
    _dim(_mesh.dimension()),
    _step(_fe_problem.timeStep()),
    _deviatroic_strain_rate(declareProperty<Real>("deviatoric_strain_rate")),
    _deviatroic_strain_rate_old(getMaterialPropertyOldByName<Real>("deviatoric_strain_rate")),
    _Cd_mat(declareProperty<Real>("Cd_mat")),
    _Cd_mat_old(getMaterialPropertyOldByName<Real>("Cd_mat")),
    _use_strain_rate_dependent_Cd(getParam<bool>("use_strain_rate_dependent_Cd")),
    _m_exponent(getParam<Real>("m_exponent")),
    _strain_rate_hat(getParam<Real>("strain_rate_hat")),
    _cd_hat(getParam<Real>("cd_hat")),
    _zero_Cd_below_threshold(getParam<bool>("zero_Cd_below_threshold")),
    //use nonlocal equivalent strain
    _use_nonlocal_eqstrain(getParam<bool>("use_nonlocal_eqstrain")),
    _eqstrain_nonlocal_old(getMaterialPropertyOldByName<Real>("eqstrain_nonlocal")),
    _nonlocal_eqstrain_blocks(getParam<std::vector<unsigned int>>("nonlocal_eqstrain_blocks")),
    //use nonlocal strain rate for Cd calculation
    _use_nonlocal_strain_rate(getParam<bool>("use_nonlocal_strain_rate")),
    _strain_rate_nonlocal_old(_use_nonlocal_strain_rate ?
        &getMaterialPropertyOld<Real>("strain_rate_nonlocal") : nullptr),
    //static solve flag
    _static_solve_flag(getParam<bool>("static_solve_flag"))
{
  // Enforce explicit block list when nonlocal eqstrain is enabled
  if (_use_nonlocal_eqstrain && _nonlocal_eqstrain_blocks.empty())
    mooseError("When 'use_nonlocal_eqstrain=true' you must provide 'nonlocal_eqstrain_blocks' (e.g., 'nonlocal_eqstrain_blocks = 100 200').");
}

void
ComputeDamageBreakageStress3DSlipWeakeningNonlocalRestart::initialSetup()
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
ComputeDamageBreakageStress3DSlipWeakeningNonlocalRestart::initQpStatefulProperties()
{
  // _elastic_strain[_qp].zero();
  // _stress[_qp].zero();
  // _deviatroic_strain_rate[_qp] = 0.0;
  // _Cd_mat[_qp] = 0.0;
}

void
ComputeDamageBreakageStress3DSlipWeakeningNonlocalRestart::computeQpStress()
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

  // std::cout << "gamma_damaged_r: " << gamma_damaged_r << std::endl;
  // std::cout << "a0: " << a0 << ", a1: " << a1 << ", a2: " << a2 << ", a3: " << a3 << std::endl;

  // if (_step == 1){
  //   setupInitial();
  //   _stress[_qp].zero();
  // }
  // else{

    /*
    compute alpha and B parameters
    */

    //compute Cd
    if (_use_strain_rate_dependent_Cd) // strain rate dependent Cd
      computeStrainRateCd();
    else // constant Cd
      _Cd_mat[_qp] = _Cd_constant;

    /* compute alpha */
    //compute forcing term
    Real alpha_forcingterm;
    Real xi_old = _xi_old[_qp];
    if (useNonlocalEqStrainHere())
    {
      xi_old = _eqstrain_nonlocal_old[_qp];
    }

    if ( xi_old >= _xi_0 && xi_old <= _xi_max ){
      alpha_forcingterm = (1 - _B_old[_qp]) * ( _Cd_mat[_qp] * _I2_old[_qp] * ( xi_old - _xi_0 ) );
    }
    else if ( xi_old < _xi_0 && xi_old >= _xi_min ){
      alpha_forcingterm = (1 - _B_old[_qp]) * ( _C1 * std::exp(_alpha_damagedvar_old[_qp]/_C2) * _I2_old[_qp] * ( xi_old - _xi_0 ) );
    }
    else{
      mooseError("xi_old is OUT-OF-RANGE!.");
    }

    //update alpha at current time
    Real alpha_out = _alpha_damagedvar_old[_qp] + _dt * alpha_forcingterm;

    //check alpha within range
    if ( alpha_out < 0 ){ alpha_out = 0.0; }
    else if ( alpha_out > 1 ){ alpha_out = 1.0; }
    else{}

    //check below initial damage (fix initial damage)
    if ( alpha_out < _initial_damage[_qp] + _damage_perturbation[_qp]){ alpha_out = _initial_damage[_qp] + _damage_perturbation[_qp]; }
    else{}

    _alpha_damagedvar[_qp] = alpha_out;

    /* compute B */
    Real C_B = _CdCb_multiplier * _Cd_mat[_qp]; //multiplier between Cd and Cb

    //compute xi_1
    Real _xi_1 = _xi_0 + sqrt( pow(_xi_0 , 2) + 2 * _shear_modulus_o / _lambda_o );

    //alphacr function
    Real alphacr;
    if ( xi_old < _xi_0 ){ alphacr = 1.0;}
    else if ( xi_old > _xi_0 && xi_old <= _xi_1 ){ alphacr = alphacr_root1(xi_old,gamma_damaged_r);}
    else if ( xi_old > _xi_1 && xi_old <= _xi_max ){ alphacr = alphacr_root2(xi_old,gamma_damaged_r); }
    else{std::cout<<"xi: "<<xi_old<<std::endl;mooseError("xi exceeds the maximum allowable range!");}

    //compute forcing func
    Real Prob = 1.0 / ( std::exp( (alphacr - _alpha_damagedvar_old[_qp]) / _beta_width ) + 1.0 );
    Real B_forcingterm;
    if ( xi_old >= _xi_d && xi_old <= _xi_max ){
      B_forcingterm = 1.0 * C_B * Prob * (1-_B_old[_qp]) * _I2_old[_qp] * (xi_old - _xi_d); //could heal if xi < xi_0
    }
    else if ( xi_old < _xi_d && xi_old >= _xi_min ){
      B_forcingterm = 1.0 * _CBH_constant * _I2_old[_qp] * ( xi_old - _xi_d ); //close healing
    }
    else{
      mooseError("xi_old is OUT-OF-RANGE!.");
    }

    Real B_out = _B_old[_qp] + _dt * B_forcingterm;

    //check breakage within range
    if ( B_out < 0 ){ B_out = 0.0; }
    else if ( B_out > 1 ){ B_out = 1.0; }
    else{}

    //check below initial damage (fix initial damage)
    if ( B_out < _initial_breakage[_qp] ){ B_out = _initial_breakage[_qp]; }
    else{}

    //save alpha and B
    _B[_qp] = B_out;

    //lambda, shear_modulus, gamma_damaged are updated
    Real lambda_out = _lambda_o;
    Real shear_modulus_out = _shear_modulus_o + alpha_out * _xi_0 * gamma_damaged_r;
    Real gamma_damaged_out = alpha_out * gamma_damaged_r;

    //save
    _lambda[_qp] = lambda_out;
    _shear_modulus[_qp] = shear_modulus_out;
    _gamma_damaged[_qp] = gamma_damaged_out;

    /* compute strain */
    RankTwoTensor eps_p = _eps_p_old[_qp] + _dt * _C_g * std::pow(_B_old[_qp],_m1) * _sigma_d_old[_qp];
    RankTwoTensor eps_t_inc = _mechanical_strain[_qp] - _mechanical_strain_old[_qp];
    RankTwoTensor eps_total = _eps_total_old[_qp] + eps_t_inc;
    RankTwoTensor eps_e = eps_total - eps_p;

    const Real epsilon = 1e-12;
    Real I1 = epsilon + eps_e(0,0) + eps_e(1,1) + eps_e(2,2);
    Real I2 = epsilon + eps_e(0,0) * eps_e(0,0) + eps_e(1,1) * eps_e(1,1) + eps_e(2,2) * eps_e(2,2) + 2 * eps_e(0,1) * eps_e(0,1) + 2 * eps_e(0,2) * eps_e(0,2) + 2 * eps_e(1,2) * eps_e(1,2);
    Real xi = I1/std::sqrt(I2);

    //Represent sigma (solid(s) + granular(b))
    RankTwoTensor sigma_s;
    RankTwoTensor sigma_b;
    RankTwoTensor sigma_total;
    RankTwoTensor sigma_d;
    const auto I = RankTwoTensor::Identity();

    /* Compute stress */
    sigma_s = (lambda_out - gamma_damaged_out / xi) * I1 * RankTwoTensor::Identity() + (2 * shear_modulus_out - gamma_damaged_out * xi) * eps_e;
    sigma_b = (2 * a2 + a1 / xi + 3 * a3 * xi) * I1 * RankTwoTensor::Identity() + (2 * a0 + a1 * xi - a3 * std::pow(xi, 3)) * eps_e;
    sigma_total = (1 - B_out) * sigma_s + B_out * sigma_b;

    sigma_d = sigma_total - 0.3333 * (sigma_total(0,0) + sigma_total(1,1) + sigma_total(2,2)) * I;

    _eps_total[_qp] = eps_p + eps_e;
    _eps_p[_qp] = eps_p;
    _eps_e[_qp] = eps_e;
    _I1[_qp] = I1;
    _I2[_qp] = I2;
    _xi[_qp] = xi;
    _sigma_d[_qp] = sigma_d;

    // Rotate the stress state to the current configuration
    // Here the stress increments are feed into the stress tensor
    _stress[_qp] = sigma_total - _sts_initial_tensor_old[_qp];

    // Also save the total stress tensor
    _sts_total[_qp] = sigma_total;

    // Always take the old value of initial stress tensor
    _sts_initial_tensor[_qp] = _sts_initial_tensor_old[_qp];

    // Assign value for elastic strain, which is equal to the mechanical strain
    _elastic_strain[_qp] = eps_e; //- _static_initial_strain_tensor[_qp];

    // Compute tangent
    RankFourTensor tangent;
    computeQpTangentModulus(tangent,I1,I2,xi,eps_e,
                            a0,a1,a2,a3,gamma_damaged_r);
    _Jacobian_mult[_qp] = tangent;

    //Compute deviatoric strain rate tensor
    computeDeviatroicStrainRateTensor();
  // }

}

Real
ComputeDamageBreakageStress3DSlipWeakeningNonlocalRestart::computegammar()
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
ComputeDamageBreakageStress3DSlipWeakeningNonlocalRestart::computecoefficients(Real gamma_damaged_r)
{

  //compute xi_1
  Real _xi_1 = _xi_0 + sqrt( pow(_xi_0 , 2) + 2 * _shear_modulus_o / _lambda_o );

  // std::cout << "xi_1: " << _xi_1 << std::endl;

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
ComputeDamageBreakageStress3DSlipWeakeningNonlocalRestart::alphacr_root1(Real xi, Real gamma_damaged_r) {
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
ComputeDamageBreakageStress3DSlipWeakeningNonlocalRestart::alphacr_root2(Real xi, Real gamma_damaged_r) {
    return 2 * _shear_modulus_o / (gamma_damaged_r * (xi - 2 * _xi_0));
}

void
ComputeDamageBreakageStress3DSlipWeakeningNonlocalRestart::computeQpTangentModulus(RankFourTensor & tangent,
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
  Real shear_modulus_out = _shear_modulus_o + _alpha_damagedvar[_qp] * _xi_0 * gamma_damaged_r;
  Real gamma_damaged_out = _alpha_damagedvar[_qp] * gamma_damaged_r;

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
  tangent = dSsdE * (1.0 - _B[_qp]) + dSbdE * _B[_qp];

}

void
ComputeDamageBreakageStress3DSlipWeakeningNonlocalRestart::setupInitial()
{

  //if no static solve is done before, then the initial strain is computed through hook's law
  //with provided initial stress
  //the initial damage and breakage should be set to zero
  if(!_static_solve_flag){
    /// lambda (first lame const)
    _lambda[_qp] = _lambda_o;
    /// mu (shear modulus)
    _shear_modulus[_qp] = _shear_modulus_o;
    /// gamma_damaged (damage modulus)
    _gamma_damaged[_qp] = 0.0;

    //allpha, B
    _alpha_damagedvar[_qp] = 0.0;
    _B[_qp] = 0.0;

    //Get stress components
    RankTwoTensor stress_initial = _static_initial_stress_tensor[_qp];
    // RankTwoTensor strain_initial = _static_initial_strain_tensor[_qp];

    //Extract stress components from stress tensor
    Real sts11_init = stress_initial(0,0);
    Real sts22_init = stress_initial(1,1);
    Real sts33_init = stress_initial(2,2);
    Real sts12_init = stress_initial(0,1);
    Real sts13_init = stress_initial(0,2);
    Real sts23_init = stress_initial(1,2);

    //Convert (lambda_o,shear_modulus_o) to (youngs_modulus_o,poisson_ratio_o)
    Real youngs_modulus_o = _shear_modulus_o * ( 3 * _lambda_o + 2 * _shear_modulus_o ) / ( _lambda_o + _shear_modulus_o );
    Real poisson_ratio_o = _lambda_o / ( 2 * ( _lambda_o + _shear_modulus_o ));

    //Compute strain components using Hooke's Law
    Real eps11_init, eps22_init, eps33_init, eps12_init, eps13_init, eps23_init;

    if (_dim == 2) {
      // 2D PLANE STRAIN: eps_zz = 0
      // For plane strain: sigma_zz = nu*(sigma_xx + sigma_yy)
      // In-plane strain-stress relations:
      //   eps_xx = (1+nu)/E * [(1-nu)*sigma_xx - nu*sigma_yy]
      //   eps_yy = (1+nu)/E * [(1-nu)*sigma_yy - nu*sigma_xx]
      //   eps_xy = (1+nu)/E * sigma_xy

      eps11_init = (1.0 + poisson_ratio_o) / youngs_modulus_o * ((1.0 - poisson_ratio_o) * sts11_init - poisson_ratio_o * sts22_init);
      eps22_init = (1.0 + poisson_ratio_o) / youngs_modulus_o * ((1.0 - poisson_ratio_o) * sts22_init - poisson_ratio_o * sts11_init);
      eps12_init = (1.0 + poisson_ratio_o) / youngs_modulus_o * sts12_init;
      eps13_init = 0.0; // plane strain: out-of-plane shear = 0
      eps23_init = 0.0; // plane strain: out-of-plane shear = 0
      eps33_init = 0.0; // plane strain: eps_zz = 0
    }
    else {
      // 3D ISOTROPIC
      eps11_init = 1.0 / youngs_modulus_o * ( sts11_init - poisson_ratio_o * ( sts22_init + sts33_init ) );
      eps22_init = 1.0 / youngs_modulus_o * ( sts22_init - poisson_ratio_o * ( sts11_init + sts33_init ) );
      eps33_init = 1.0 / youngs_modulus_o * ( sts33_init - poisson_ratio_o * ( sts11_init + sts22_init ) );
      eps12_init = 1.0 / youngs_modulus_o * ( ( 1 + poisson_ratio_o ) * sts12_init );
      eps13_init = 1.0 / youngs_modulus_o * ( ( 1 + poisson_ratio_o ) * sts13_init );
      eps23_init = 1.0 / youngs_modulus_o * ( ( 1 + poisson_ratio_o ) * sts23_init );
    }

    //Compute xi, I1, I2
    Real I1_init = eps11_init + eps22_init + eps33_init;
    Real I2_init = eps11_init * eps11_init + eps22_init * eps22_init + eps33_init * eps33_init + 2 * eps12_init * eps12_init + 2 * eps13_init * eps13_init + 2 * eps23_init * eps23_init;
    Real xi_init = I1_init / sqrt( I2_init );

    //Compute eps
    //eps_p
    _eps_p[_qp](0,0) = 0.0; _eps_p[_qp](0,1) = 0.0; _eps_p[_qp](0,2) = 0.0;
    _eps_p[_qp](1,0) = 0.0; _eps_p[_qp](1,1) = 0.0; _eps_p[_qp](1,2) = 0.0;
    _eps_p[_qp](2,0) = 0.0; _eps_p[_qp](2,1) = 0.0; _eps_p[_qp](2,2) = 0.0;
    //eps_e
    _eps_e[_qp](0,0) = eps11_init; _eps_e[_qp](0,1) = eps12_init; _eps_e[_qp](0,2) = eps13_init;
    _eps_e[_qp](1,0) = eps12_init; _eps_e[_qp](1,1) = eps22_init; _eps_e[_qp](1,2) = eps23_init;
    _eps_e[_qp](2,0) = eps13_init; _eps_e[_qp](2,1) = eps23_init; _eps_e[_qp](2,2) = eps33_init;
    //eps_total
    _eps_total[_qp](0,0) = eps11_init; _eps_total[_qp](0,1) = eps12_init; _eps_total[_qp](0,2) = eps13_init;
    _eps_total[_qp](1,0) = eps12_init; _eps_total[_qp](1,1) = eps22_init; _eps_total[_qp](1,2) = eps23_init;
    _eps_total[_qp](2,0) = eps13_init; _eps_total[_qp](2,1) = eps23_init; _eps_total[_qp](2,2) = eps33_init;
    //sts_total
    _sts_total[_qp] = stress_initial;

    //I1
    _I1[_qp] = I1_init;
    //I2
    _I2[_qp] = I2_init;
    //xi
    _xi[_qp] = xi_init;

    //Compute sigma_d (deviatoric stress)
    const auto I = RankTwoTensor::Identity();
    RankTwoTensor sigma_d = stress_initial - 0.3333 * (stress_initial(0,0) + stress_initial(1,1) + stress_initial(2,2)) * I;
    _sigma_d[_qp] = sigma_d;

    // Also save in the sts_initial_tensor
    _sts_initial_tensor[_qp] = stress_initial;

    // Assign value for elastic strain, which is equal to the mechanical strain
    _elastic_strain[_qp] = _eps_e[_qp]; //- _static_initial_strain_tensor[_qp];
  }
  //if the static solve is performed, then the strain should read initial strain tensor
  // and the damage and breakage could be non-zero and the modulus should be updated as well
  else{
    Real gamma_damaged_r = computegammar();
    std::vector<Real> avec = computecoefficients(gamma_damaged_r);
    Real a0 = avec[0];
    Real a1 = avec[1];
    Real a2 = avec[2];
    Real a3 = avec[3];

    //alpha, B (initialize from initial damage and breakage)
    _alpha_damagedvar[_qp] = _initial_damage[_qp];
    _B[_qp] = _initial_breakage[_qp];

    /// lambda (first lame const)
    _lambda[_qp] = _lambda_o;
    /// mu (shear modulus)
    _shear_modulus[_qp] = _shear_modulus_o + _initial_damage[_qp] * _xi_0 * gamma_damaged_r;
    /// gamma_damaged (damage modulus)
    _gamma_damaged[_qp] = _initial_damage[_qp] * gamma_damaged_r;

    RankTwoTensor eps_e = _static_initial_strain_tensor[_qp];

    const Real epsilon = 1e-12;
    Real I1 = epsilon + eps_e(0,0) + eps_e(1,1) + eps_e(2,2);
    Real I2 = epsilon + eps_e(0,0) * eps_e(0,0) + eps_e(1,1) * eps_e(1,1) + eps_e(2,2) * eps_e(2,2) + 2 * eps_e(0,1) * eps_e(0,1) + 2 * eps_e(0,2) * eps_e(0,2) + 2 * eps_e(1,2) * eps_e(1,2);
    Real xi = I1/std::sqrt(I2);

    //Represent sigma (solid(s) + granular(b))
    RankTwoTensor sigma_s;
    RankTwoTensor sigma_b;
    RankTwoTensor sigma_total;
    RankTwoTensor sigma_d;
    const auto I = RankTwoTensor::Identity();

    /* Compute stress */
    sigma_s = (_lambda[_qp] - _gamma_damaged[_qp] / xi) * I1 * RankTwoTensor::Identity() + (2 * _shear_modulus[_qp] - _gamma_damaged[_qp] * xi) * eps_e;
    sigma_b = (2 * a2 + a1 / xi + 3 * a3 * xi) * I1 * RankTwoTensor::Identity() + (2 * a0 + a1 * xi - a3 * std::pow(xi, 3)) * eps_e;
    sigma_total = (1 - _B[_qp]) * sigma_s + _B[_qp] * sigma_b;

    sigma_d = sigma_total - 0.3333 * (sigma_total(0,0) + sigma_total(1,1) + sigma_total(2,2)) * I;

    _eps_total[_qp] = eps_e;
    _eps_p[_qp].zero(); // Initialize plastic strain to zero
    _eps_e[_qp] = eps_e;
    _I1[_qp] = I1;
    _I2[_qp] = I2;
    _xi[_qp] = xi;
    _sigma_d[_qp] = sigma_d;

    // Rotate the stress state to the current configuration
    // Here the stress increments are feed into the stress tensor
    //_stress[_qp] = sigma_total - _static_initial_stress_tensor[_qp];

    // Also save the total stress tensor
    _sts_total[_qp] = sigma_total;

    // Also save in the sts_initial_tensor
    _sts_initial_tensor[_qp] = sigma_total;

    // Assign value for elastic strain, which is equal to the mechanical strain
    _elastic_strain[_qp] = eps_e; //- _static_initial_strain_tensor[_qp];
  }
}

void
ComputeDamageBreakageStress3DSlipWeakeningNonlocalRestart::computeDeviatroicStrainRateTensor()
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
ComputeDamageBreakageStress3DSlipWeakeningNonlocalRestart::computeStrainRateCd()
{
  //_m_exponent: constant value - default value = 0.8
  //_strain_rate_hat: constant value - default value = 1e-4
  //_cd_hat: constant value - default value = 1
  //_strain_rate: deviatoric strain rate, variable value passed from main app

  // Use NONLOCAL strain rate if enabled, otherwise use LOCAL strain rate
  Real effective_strain_rate;

  if (_use_nonlocal_strain_rate && useNonlocalEqStrainHere())
  {
    // Use nonlocal averaged strain rate (mesh-independent)
    effective_strain_rate = (*_strain_rate_nonlocal_old)[_qp];
  }
  else
  {
    // Use local strain rate (mesh-dependent due to dt)
    effective_strain_rate = _deviatroic_strain_rate_old[_qp];
  }

  // Compute Cd based on effective strain rate
  if (effective_strain_rate < _strain_rate_hat){
    // if deviatoric strain rate is less than strain_rate_hat, Cd = 0 (optional) or Cd_hat (default)
    _Cd_mat[_qp] = _zero_Cd_below_threshold ? 0.0 : _cd_hat;
  }
  else{
    _Cd_mat[_qp] = pow(10, 1 + _m_exponent * std::log10(effective_strain_rate/_strain_rate_hat)) * _cd_hat;
  }
}
