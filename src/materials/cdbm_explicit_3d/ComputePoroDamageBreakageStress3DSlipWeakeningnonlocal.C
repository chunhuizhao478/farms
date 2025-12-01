//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ComputePoroDamageBreakageStress3DSlipWeakeningnonlocal.h"
#include "NestedSolve.h"
#include "FEProblem.h"

registerMooseObject("farmsApp", ComputePoroDamageBreakageStress3DSlipWeakeningnonlocal);

InputParameters
ComputePoroDamageBreakageStress3DSlipWeakeningnonlocal::validParams()
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

  //Poroelastic properties
  params.addParam<Real>("permeability_solid_o", "permeability of solid meterial");
  params.addParam<Real>("initial_viscosity_fluid", "fluid viscosity");
  params.addParam<Real>("solid_bulk_modulus_s", "solid bulk modulus of solid grains");
  params.addParam<Real>("solid_bulk_modulus_g", "solid bulk modulus of granular material");
  params.addParam<Real>("fluid_bulk_modulus", "fluid bulk modulus"); 
  params.addParam<Real>("porosity_solid_o", "initial prosoity of solid phase"); 
  params.addParam<Real>("permeability_evolution_with_damage", "parameter for permeability evolution with damage"); 
  params.addParam<Real>("initial_grain_size", "initial harmonic mean grain size"); 
  params.addParam<Real>("ultimate_grain_size", "ultimate harmonic mean grain size");
  
  //Porepressure variable
  params.addCoupledVar("porepressure", 0.0, "The pore pressure variable");
  
  //Dilatancy parameters
  params.addParam<bool>("use_dilatancy", false, "Flag to use dilatancy variable evolution");
  params.addParam<Real>("anand_param_go_mat",0,"Dilatancy parameter go");
  params.addParam<Real>("anand_param_eta_cv_mat",0,"Dilatancy parameter eta_cv");
  params.addParam<Real>("anand_param_p_mat",0,"Dilatancy parameter p");

  //strain rate dependent Cd parameters
  params.addParam<bool>("use_strain_rate_dependent_Cd", false,
                        "Use strain rate dependent Cd (default: false)");
  params.addParam<Real>( "m_exponent", 0.8, "strain rate dependent parameters");
  params.addParam<Real>( "strain_rate_hat", 1e-4, "strain rate dependent parameters");
  params.addParam<Real>( "cd_hat", 1.0, "strain rate dependent parameters");
  params.addParam<bool>("zero_Cd_below_threshold", false,
                        "If true, set Cd = 0 when deviatoric strain rate < strain_rate_hat; otherwise use cd_hat.");

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

ComputePoroDamageBreakageStress3DSlipWeakeningnonlocal::ComputePoroDamageBreakageStress3DSlipWeakeningnonlocal(const InputParameters & parameters)
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
    _initial_porepressure(getMaterialProperty<Real>("initial_porepressure")),
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
    _Biot_coeff_s(declareProperty<Real>("Biot_coefficient_solid")),
    _Biot_coeff_g(declareProperty<Real>("Biot_coefficient_granular")),
    _Biot_modulus_s(declareProperty<Real>("Biot_modulus_solid")),
    _Biot_modulus_g(declareProperty<Real>("Biot_modulus_granular")),
    _biot_coeff_eff(declareProperty<Real>("biot_coefficient_effective")),
    _Biot_modulus_eff(declareProperty<Real>("Biot_modulus_effective")),
    _fluid_solid_coupling(declareProperty<Real>("fluid_solid_coupling")),
    _perm_s(declareProperty<Real>("permeability_solid")),
    _perm_g(declareProperty<Real>("permeability_granular")),
    _perm_cr(declareProperty<Real>("permeability_critical")),
    _phi_cr(declareProperty<Real>("porosity_critical")),
    _phi_p(declareProperty<Real>("plastic_porosity")),
    _permeability_solid_o(getParam<Real>("permeability_solid_o")),
    _solid_bulk_modulus_s(getParam<Real>("solid_bulk_modulus_s")),
    _solid_bulk_modulus_g(getParam<Real>("solid_bulk_modulus_g")),
    _fluid_bulk_modulus(getParam<Real>("fluid_bulk_modulus")),
    _porosity_solid_o(getParam<Real>("porosity_solid_o")),
    _initial_viscosity_fluid(getParam<Real>("initial_viscosity_fluid")),
    _b(getParam<Real>("permeability_evolution_with_damage")),
    _DHo(getParam<Real>("initial_grain_size")),
    _DHu(getParam<Real>("ultimate_grain_size")),
    _pore_pressure(coupledValue("porepressure")),
    _pore_pressure_old(coupledValueOld("porepressure")),
    _stress_off_diag_jacobian(declareProperty<RankTwoTensor>(_base_name + "stress_off_diag_jacobian")),
    // Add option to add dilatancy/compaction effect //Follow paper Section 7.1
    _use_dilatancy(getParam<bool>("use_dilatancy")),
    _eta(declareProperty<Real>(_base_name + "plastic_volume_change")),
    _eta_old(getMaterialPropertyOldByName<Real>(_base_name + "plastic_volume_change")),
    _dilatancy_function_beta(declareProperty<Real>(_base_name + "dilatancy_function_beta")),
    _shear_rate_nu(declareProperty<RankTwoTensor>(_base_name + "shear_rate_nu")),
    _anand_param_go_mat(getParam<Real>("anand_param_go_mat")),
    _anand_param_eta_cv_mat(getParam<Real>("anand_param_eta_cv_mat")),
    _anand_param_p_mat(getParam<Real>("anand_param_p_mat")),
    _deps_p_dp(declareProperty<RankTwoTensor>("deps_p_dp")),
    _deps_p_deps(declareProperty<RankFourTensor>("deps_p_deps")),
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
ComputePoroDamageBreakageStress3DSlipWeakeningnonlocal::initialSetup()
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
ComputePoroDamageBreakageStress3DSlipWeakeningnonlocal::initQpStatefulProperties()
{
  _elastic_strain[_qp].zero();
  _stress[_qp].zero();
  _I1[_qp] = 0.0;
  _deviatroic_strain_rate[_qp] = 0.0;
  _Cd_mat[_qp] = 0.0;
  _eta[_qp] = 0.0;
  _dilatancy_function_beta[_qp] = 0.0;

}

void
ComputePoroDamageBreakageStress3DSlipWeakeningnonlocal::computeQpStress()
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

  if (_step == 1){
    setupInitial();
    _stress[_qp].zero();
    _I1[_qp] = 0.0;
  }
  else{
    
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

    // Get old deviatoric stress
    RankTwoTensor sigma_d_old_tensor = _sigma_d_old[_qp];

    // Get the norm of deviatoric stress scalar
    Real sigma_d_norm = 0.0;
    for (unsigned int p = 0; p < 3; p++){
        for (unsigned int q = 0; q < 3; q++){
            sigma_d_norm += sigma_d_old_tensor(p,q) * sigma_d_old_tensor(p,q);
        }
    }
    sigma_d_norm = std::sqrt(sigma_d_norm);

    // Get deviatoric stress direction
    RankTwoTensor N; 
    N.zero();

    // Epsilon to avoid division by zero
    if (sigma_d_norm != 0.0){
        // Compute deviatoric stress direction
        for (unsigned int p = 0; p < 3; p++){
            for (unsigned int q = 0; q < 3; q++){
                N(p,q) = sigma_d_old_tensor(p,q) / sigma_d_norm;
            }
        }
    }

    // Define equivalent plastic strain rate
    RankTwoTensor eps_p_dot;

    if (_use_dilatancy)
    {
        _shear_rate_nu[_qp] = _C_g * std::pow(_B_old[_qp], _m1) * std::pow(sigma_d_norm, _m2) * N;
        
        _eta[_qp] = _eta_old[_qp] + _dilatancy_function_beta[_qp] * _C_g * std::pow(_B_old[_qp], _m1) * std::pow(sigma_d_norm, _m2) * _dt;
        
        _dilatancy_function_beta[_qp] = _anand_param_go_mat * std::pow(1 - _eta[_qp] / _anand_param_eta_cv_mat, _anand_param_p_mat );  
        
        // Plastic strain rate with dilatancy (deviatoric + volumetric parts)
        eps_p_dot = _shear_rate_nu[_qp] + _C_g * std::pow(_B_old[_qp], _m1) * std::pow(sigma_d_norm, _m2) * _dilatancy_function_beta[_qp]/ 3.0 * RankTwoTensor::Identity();
    }
    else
    {
       _shear_rate_nu[_qp] = _C_g * std::pow(_B_old[_qp], _m1) * std::pow(sigma_d_norm, _m2) * N;
        
        // Plastic strain rate without dilatancy (deviatoric only)
        eps_p_dot = _shear_rate_nu[_qp];
    }

    // Update plastic strain using backward Euler (small strain)
    RankTwoTensor eps_p = _eps_p_old[_qp] + _dt * eps_p_dot;

    RankTwoTensor eps_t_inc = _mechanical_strain[_qp] - _mechanical_strain_old[_qp];
    RankTwoTensor eps_total = _eps_total_old[_qp] + eps_t_inc;
    RankTwoTensor eps_e = eps_total - eps_p;

    const Real epsilon = 1e-12;
    Real I1 = epsilon + eps_e(0,0) + eps_e(1,1) + eps_e(2,2);
    Real I2 = epsilon + eps_e(0,0) * eps_e(0,0) + eps_e(1,1) * eps_e(1,1) + eps_e(2,2) * eps_e(2,2) + 2 * eps_e(0,1) * eps_e(0,1) + 2 * eps_e(0,2) * eps_e(0,2) + 2 * eps_e(1,2) * eps_e(1,2);
    Real xi = I1/std::sqrt(I2);

    /* poroelastic properties solid and granular phase */

    // Solid bulk modulus (constant for solid grains)
    Real K_s = _solid_bulk_modulus_s;

    // Solid bulk modulus after crushing
    Real K_s_crushed = _solid_bulk_modulus_g;
    
    // Fluid bulk modulus
    Real K_f = _fluid_bulk_modulus;
    
    // Compute drained bulk modulus K_d of solid phase
    Real K_d = _lambda[_qp] + (2.0/3.0) * _shear_modulus[_qp] - (2.0/3.0) * _gamma_damaged[_qp] * xi;

    // Compute bulk modulus for granular phase
    Real K_d_granular = 2 * a2 + a3 * (6.0 - (4.0/3.0) *  xi * xi ) * xi 
                        + (2.0/3.0) * a0 + (2.0/3.0) * a1 * xi;

    // Solid bulk modulus for granular material 
    Real K_s_granular = (1 - _B[_qp]) * K_s + _B[_qp] * K_s_crushed;

    // Compute Biot coefficient for solid phase
    Real alpha_s = 1.0 - K_d/K_s;

    // Compute Biot coefficient for granular phase
    Real alpha_g = 1.0 - K_d_granular/K_s_granular;
    
    // Compute porosity evolution for solid phase
    // Real porosity_s = 1 - (1 - _porosity_solid_o) * exp(-I1);
     Real porosity_s = _porosity_solid_o;

    // Compute plastic porosity evolution 
    Real porosity_p = eps_p(0,0) + eps_p(1,1) + eps_p(2,2);

    // Compute total porosity
    Real porosity = porosity_s + porosity_p;

    // Compute Biot modulus for solid phase
    Real one_over_Storage_s = (K_s*K_f)/(porosity_s * K_f + (alpha_s - porosity_s) * K_s);
    
    // Compute Biot modulus for solid phase
    Real one_over_Storage_g = (K_s_granular*K_f)/(porosity * K_f + (alpha_g - porosity) * K_s_granular);

    // Compute permeability for solid phase
    Real perm_s = _permeability_solid_o * pow(porosity_s/_porosity_solid_o, 3.0) * exp(_b * alpha_out);
    
    // Determine critical porosity and permeability based on damage state
    const Real tolerance = 1e-12;  // Small tolerance for floating point comparison

    Real alpha_cr = alphacr_root1(xi, gamma_damaged_r);

    if (std::abs(_alpha_damagedvar[_qp] - alpha_cr) < tolerance) {
      // At critical damage: use current solid phase properties
      _phi_cr[_qp] = porosity_s;
      _perm_cr[_qp] = perm_s;
    } else {
      // Below critical damage: use initial solid properties
      _phi_cr[_qp] = _porosity_solid_o;
      _perm_cr[_qp] = _permeability_solid_o;
    }

    // Compute harmonic mean grain size for current distribution
    Real DH = (1 - _B[_qp]) * _DHo + _B[_qp] * _DHu;
   
    // Compute permeability for granular phase
    Real perm_g = _perm_cr[_qp] * pow(porosity/_phi_cr[_qp], 3.0) * pow(DH/_DHo, 2.0);
    // Real perm_g = _perm_cr[_qp] * pow(porosity/_phi_cr[_qp], 3.0);

    // Save solid phase properties
    _Biot_coeff_s[_qp] = alpha_s;
    _Biot_modulus_s[_qp] = one_over_Storage_s;
    _perm_s[_qp] = (1 - _B[_qp]) * perm_s / _initial_viscosity_fluid;

    // Save granular phase properties
    _Biot_coeff_g[_qp] = alpha_g;
    _Biot_modulus_g[_qp] = one_over_Storage_g;
    _phi_p[_qp] = porosity_p;
    _perm_g[_qp] = _B[_qp] * perm_g /_initial_viscosity_fluid;
    
    //Represent sigma (solid(s) + granular(b))
    RankTwoTensor sigma_s;
    RankTwoTensor sigma_b;
    RankTwoTensor fluid_contribution;
    RankTwoTensor sigma_total;
    RankTwoTensor sigma_eff;
    RankTwoTensor sigma_d_eff;

    Real term11 = (1 - _B[_qp]) * _B[_qp] * _Biot_modulus_s[_qp] * _Biot_modulus_g[_qp] * std::pow((_Biot_coeff_s[_qp] - _Biot_coeff_g[_qp]), 2);
    Real term22 = (1 - _B[_qp]) * _Biot_coeff_s[_qp] * _Biot_modulus_s[_qp] + _B[_qp] * _Biot_coeff_g[_qp] * _Biot_modulus_g[_qp];
    Real term33 = (1 - _B[_qp]) * _Biot_modulus_s[_qp] + _B[_qp] * _Biot_modulus_g[_qp];

    _Biot_modulus_eff[_qp] = term33;
    _fluid_solid_coupling[_qp] = term22 / term33;
    _biot_coeff_eff[_qp] = term22 / term33;
  
    const auto I = RankTwoTensor::Identity();

    /* Compute stress */
    sigma_s = (lambda_out - gamma_damaged_out / xi) * I1 * RankTwoTensor::Identity() + (2 * shear_modulus_out - gamma_damaged_out * xi) * eps_e;
    sigma_b = (2 * a2 + a1 / xi + 3 * a3 * xi) * I1 * RankTwoTensor::Identity() + (2 * a0 + a1 * xi - a3 * std::pow(xi, 3)) * eps_e;
    fluid_contribution = term11 / term33 * I1 * RankTwoTensor::Identity() - term22 / term33 * (_initial_porepressure[_qp] + _pore_pressure[_qp])  * RankTwoTensor::Identity();
    sigma_total = (1 - B_out) * sigma_s + B_out * sigma_b + fluid_contribution;

    sigma_eff = sigma_total +  (_initial_porepressure[_qp] + _pore_pressure[_qp]) * RankTwoTensor::Identity();
    
    sigma_d_eff = sigma_eff - 0.3333 * (sigma_eff(0,0) + sigma_eff(1,1) + sigma_eff(2,2)) * I;

    _eps_total[_qp] = eps_p + eps_e;
    _eps_p[_qp] = eps_p;
    _eps_e[_qp] = eps_e;
    _I1[_qp] = eps_t_inc(0,0)+eps_t_inc(1,1)+eps_t_inc(2,2);
    _I2[_qp] = I2;
    _xi[_qp] = xi;
    _sigma_d[_qp] = sigma_d_eff;

    // Rotate the stress state to the current configuration
    // Here the stress increments are feed into the stress tensor
    _stress[_qp] = sigma_total - _sts_initial_tensor_old[_qp];

    // Also save the total stress tensor
    _sts_total[_qp] = sigma_total;

    // Always take the old value of initial stress tensor
    _sts_initial_tensor[_qp] = _sts_initial_tensor_old[_qp];

    // Assign value for elastic strain, which is equal to the mechanical strain
    _elastic_strain[_qp] = eps_e ; //- _static_initial_strain_tensor[_qp];

    // Compute tangent
    RankFourTensor tangent;
    computeQpTangentModulus(tangent,I1,I2,xi,eps_e,
                            a0,a1,a2,a3,gamma_damaged_r);
    _Jacobian_mult[_qp] = tangent;

    //Compute deviatoric strain rate tensor
    computeDeviatroicStrainRateTensor();

    // Compute simplified derivatives using finite differences
    computeSimplifiedPlasticDerivatives(eps_p, eps_p_dot);
  }

}

Real 
ComputePoroDamageBreakageStress3DSlipWeakeningnonlocal::computegammar()
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
ComputePoroDamageBreakageStress3DSlipWeakeningnonlocal::computecoefficients(Real gamma_damaged_r)
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
ComputePoroDamageBreakageStress3DSlipWeakeningnonlocal::alphacr_root1(Real xi, Real gamma_damaged_r) {
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
ComputePoroDamageBreakageStress3DSlipWeakeningnonlocal::alphacr_root2(Real xi, Real gamma_damaged_r) {
    return 2 * _shear_modulus_o / (gamma_damaged_r * (xi - 2 * _xi_0));
}

void
ComputePoroDamageBreakageStress3DSlipWeakeningnonlocal::computeQpTangentModulus(RankFourTensor & tangent, 
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

  Real term11 = (1 - _B[_qp]) * _B[_qp] * _Biot_modulus_s[_qp] * _Biot_modulus_g[_qp] * std::pow((_Biot_coeff_s[_qp] - _Biot_coeff_g[_qp]), 2);
  Real term22 = (1 - _B[_qp]) * _Biot_coeff_s[_qp] * _Biot_modulus_s[_qp] + _B[_qp] * _Biot_coeff_g[_qp] * _Biot_modulus_g[_qp];
  Real term33 = (1 - _B[_qp]) * _Biot_modulus_s[_qp] + _B[_qp] * _Biot_modulus_g[_qp];
  

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

  // Compute fluid tangent contribution: d(fluid_contribution)/dE
  // fluid_contribution = term11/term33 * I1 * Identity - term22/term33 * p * Identity
  // d(fluid_contribution)/dE = term11/term33 * dI1/dE * Identity
  // where dI1/dE = Identity (trace operation)
    
  Real fluid_contrib = term11 / term33;
    
  RankFourTensor fluid_tangent;
  fluid_tangent.zero();
    
  // d(fluid_contribution_ij)/dE_kl = fluid_coeff * d(I1)/dE_kl * δ_ij
  // where d(I1)/dE_kl = δ_kl
  for (unsigned int i = 0; i < 3; ++i) {
    for (unsigned int j = 0; j < 3; ++j) {
      for (unsigned int k = 0; k < 3; ++k) {
        for (unsigned int l = 0; l < 3; ++l) {
            fluid_tangent(i, j, k, l) = fluid_contrib * identity(k, l) * identity(i, j);
        }
      }
    }
  }
    
  // Add fluid tangent to total tangent
  tangent += fluid_tangent;

  // Final derivative: dS/dp = -α_eff * I
  _stress_off_diag_jacobian[_qp] = -(term22 / term33) * RankTwoTensor::Identity();

}

void
ComputePoroDamageBreakageStress3DSlipWeakeningnonlocal::setupInitial()
{


  Real gamma_damaged_r = computegammar();
  std::vector<Real> avec = computecoefficients(gamma_damaged_r);
  Real a0 = avec[0];
  Real a1 = avec[1];
  Real a2 = avec[2];
  Real a3 = avec[3];

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


  /* poroelastic properties solid and granular phase */

  // Solid bulk modulus (constant for solid grains)
  Real K_s = _solid_bulk_modulus_s;

  // Fluid bulk modulus
  Real K_f = _fluid_bulk_modulus;
    
  // Compute drained bulk modulus K_d of solid phase
  Real K_d = _lambda[_qp] + (2.0/3.0) * _shear_modulus[_qp] - (2.0/3.0) * _gamma_damaged[_qp] * xi;

  // Compute Biot coefficient for solid phase
  Real alpha_s = 1.0 - K_d/K_s;
    
  // Compute porosity evolution for solid phase
  Real porosity_s = _porosity_solid_o;

  // Compute Biot modulus for solid phase
  Real one_over_Storage_s = (K_s*K_f)/(porosity_s * K_f + (alpha_s - porosity_s) * K_s);

  // Compute permeability for solid phase
  Real perm_s = _permeability_solid_o;

  // Save solid phase properties
  _Biot_coeff_s[_qp] = alpha_s;
  _Biot_modulus_s[_qp] = one_over_Storage_s;
  _perm_s[_qp] = (1 - _B[_qp]) * perm_s / _initial_viscosity_fluid;

  Real term22 = (1 - _B[_qp]) * _Biot_coeff_s[_qp] * _Biot_modulus_s[_qp];
  Real term33 = (1 - _B[_qp]) * _Biot_modulus_s[_qp];
  
  //Represent sigma (solid(s) + granular(b))
  RankTwoTensor sigma_s;
  RankTwoTensor sigma_b;
  RankTwoTensor fluid_contribution;
  RankTwoTensor sigma_total;
  RankTwoTensor sigma_eff;
  RankTwoTensor sigma_d_eff;
  const auto I = RankTwoTensor::Identity();

  /* Compute stress */
  sigma_s = (_lambda[_qp] - _gamma_damaged[_qp] / xi) * I1 * RankTwoTensor::Identity() + (2 * _shear_modulus[_qp] - _gamma_damaged[_qp] * xi) * eps_e;
  sigma_b = (2 * a2 + a1 / xi + 3 * a3 * xi) * I1 * RankTwoTensor::Identity() + (2 * a0 + a1 * xi - a3 * std::pow(xi, 3)) * eps_e;
  fluid_contribution = - term22 / term33 * (_initial_porepressure[_qp] +_pore_pressure[_qp])  * RankTwoTensor::Identity();
  sigma_total = (1 - _B[_qp]) * sigma_s + _B[_qp] * sigma_b + fluid_contribution;

  sigma_eff = sigma_total +  _initial_porepressure[_qp] * RankTwoTensor::Identity();
    
  sigma_d_eff = sigma_eff - 0.3333 * (sigma_eff(0,0) + sigma_eff(1,1) + sigma_eff(2,2)) * I;

  _eps_total[_qp] = eps_e;
  _eps_p[_qp].zero(); // Initialize plastic strain to zero
  _eps_e[_qp] = eps_e;
  _I1[_qp] = 0;
  _I2[_qp] = I2;
  _xi[_qp] = xi;
  _sigma_d[_qp] = sigma_d_eff;

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

void
ComputePoroDamageBreakageStress3DSlipWeakeningnonlocal::computeDeviatroicStrainRateTensor()
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
ComputePoroDamageBreakageStress3DSlipWeakeningnonlocal::computeStrainRateCd()
{
  //_m_exponent: constant value - default value = 0.8
  //_strain_rate_hat: constant value - default value = 1e-4
  //_cd_hat: constant value - default value = 1
  //_strain_rate: deviatoric strain rate, variable value passed from main app
  if (_deviatroic_strain_rate_old[_qp] < _strain_rate_hat){
    // if deviatoric strain rate is less than strain_rate_hat, Cd = 0 (optional) or Cd_hat (default)
    _Cd_mat[_qp] = _zero_Cd_below_threshold ? 0.0 : _cd_hat;
  }
  else{
    _Cd_mat[_qp] = pow(10, 1 + _m_exponent * std::log10(_deviatroic_strain_rate_old[_qp]/_strain_rate_hat)) * _cd_hat;
  }
}

void
ComputePoroDamageBreakageStress3DSlipWeakeningnonlocal::computeSimplifiedPlasticDerivatives(
    const RankTwoTensor & eps_p,
    const RankTwoTensor & eps_p_dot)
{
  // Approximate ∂εᵖ/∂p using finite differences
  // Similar to the finite strain approach: dFp/dp ≈ Fp_dot / dp
  
  Real dp = _pore_pressure[_qp] - _pore_pressure_old[_qp];
  
  _deps_p_dp[_qp].zero();
  
  if (_dt != 0.0 && std::abs(dp) > 1e-12)
  {
    // Approximate: ∂εᵖ_ij/∂p ≈ ε̇ᵖ_ij / (dp/dt) = (εᵖ - εᵖ_old) / dp
    for (unsigned int i = 0; i < 3; ++i)
    {
      for (unsigned int j = 0; j < 3; ++j)
      {
        _deps_p_dp[_qp](i, j) = (eps_p(i, j) - _eps_p_old[_qp](i, j)) / dp;
      }
    }
  }
  
  // Approximate ∂εᵖ/∂ε using finite differences
  // This is more complex - we approximate based on strain increment
  
  RankTwoTensor eps_inc = _mechanical_strain[_qp] - _mechanical_strain_old[_qp];
  Real eps_inc_norm = eps_inc.L2norm();
  
  _deps_p_deps[_qp].zero();
  
  if (_dt != 0.0 && eps_inc_norm > 1e-12)
  {
    // Approximate: ∂εᵖ_ij/∂ε_kl ≈ Δεᵖ_ij / Δε_kl
    for (unsigned int i = 0; i < 3; ++i)
    {
      for (unsigned int j = 0; j < 3; ++j)
      {
        for (unsigned int k = 0; k < 3; ++k)
        {
          for (unsigned int l = 0; l < 3; ++l)
          {
            if (std::abs(eps_inc(k, l)) > 1e-12)
            {
              _deps_p_deps[_qp](i, j, k, l) = 
                  (eps_p(i, j) - _eps_p_old[_qp](i, j)) / eps_inc(k, l);
            }
          }
        }
      }
    }
  }
}