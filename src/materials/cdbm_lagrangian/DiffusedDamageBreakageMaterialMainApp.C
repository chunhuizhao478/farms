//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DiffusedDamageBreakageMaterialMainApp.h"

/**
 *  Material used in damage-breakage large deformation formulation, consider full damage evolution equation with diffusion
 *  Created by Chunhui Zhao, Dec 24th, 2024
 */
registerMooseObject("farmsApp", DiffusedDamageBreakageMaterialMainApp);

InputParameters
DiffusedDamageBreakageMaterialMainApp::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Material used in three field poro dynamics simulations");
  //input parameters
  params.addParam<Real>(                   "lambda_o", "initial lambda constant value");
  params.addParam<Real>(            "shear_modulus_o", "initial shear modulus value");
  params.addParam<Real>(                       "xi_0", "strain invariants ratio: onset of damage evolution");
  params.addParam<Real>(                       "xi_d", "strain invariants ratio: onset of breakage healing");
  params.addParam<Real>(                        "chi", "coefficient of energy ratio Fb/Fs = chi < 1");
  params.addParam<Real>(                        "C_g", "compliance or fluidity of the fine grain granular material");
  params.addParam<Real>(                         "m1", "coefficient of power law indexes");
  params.addParam<Real>(                         "m2", "coefficient of power law indexes");
  //input for poroelastic material evolution  
  params.addParam<Real>(         "permeability_solid_o", "permeability of solid meterial");
  params.addParam<Real>(         "initial_viscosity_fluid", "fluid viscosity");
  params.addParam<Real>(         "solid_bulk_modulus_s", "solid bulk modulus of solid grains");
  params.addParam<Real>(         "solid_bulk_modulus_g", "solid bulk modulus of granular material");
  params.addParam<Real>(         "fluid_bulk_modulus", "fluid bulk modulus"); 
  params.addParam<Real>(         "porosity_solid_o", "initial prosoity of solid phase"); 
  params.addParam<Real>(         "permeability_evolution_with_damage", "parameter for permeability evolution with damage"); 
  params.addParam<Real>(         "initial_grain_size", "initial harmonic mean grain size"); 
  params.addParam<Real>(         "ultimate_grain_size", "ultimate harmonic mean grain size"); 
  //input coupled variables from main app
  params.addCoupledVar("structural_stress_coefficient", "structral_stress_coefficient");
  params.addCoupledVar("alpha_damagedvar_aux", "second_elastic_strain_invariant");
  params.addCoupledVar("B_damagedvar_aux", "strain_invariant_ratio");
  //build L matrix
  params.addCoupledVar(           "vel_x", "velocity in x direction"); //to build L matrix
  params.addCoupledVar(           "vel_y", "velocity in y direction"); //to build L matrix
  params.addCoupledVar(           "vel_z", "velocity in z direction"); //to build L matrix  
  //use spatial cg
  params.addParam<bool>("use_spatial_cg", false, "use spatial cg");
  params.addCoupledVar("cg_aux", "cg_aux");
  //use state dependent variables
  params.addParam<bool>("use_state_var_evolution", false, "Flag to use state variable evolution");
  params.addParam<Real>("const_A", -1.0,"Constant A value, A = a * sigma_N");
  params.addParam<Real>("const_B", -1.0,"Constant B value, B = b * sigma_N");
  params.addParam<Real>("const_theta_o", -1.0,"Constant theta_o value");
  params.addParam<Real>("initial_theta0", -1.0,"Initial theta0 value");
  params.addParam<std::string>("base_name", "", "Optional parameter for multiple material systems");
  return params;
}
DiffusedDamageBreakageMaterialMainApp::DiffusedDamageBreakageMaterialMainApp(const InputParameters & parameters)
  : Material(parameters),
  //declare properties
  //--------------------------------------------------------------//
  _gamma_damaged_r(declareProperty<Real>("gamma_damaged_r")),
  _alpha_damagedvar(declareProperty<Real>("alpha_damagedvar")),
  _B_damagedvar(declareProperty<Real>("B_damagedvar")),
  _lambda(declareProperty<Real>("lambda_const")),
  _shear_modulus(declareProperty<Real>("shear_modulus")),
  _damaged_modulus(declareProperty<Real>("damaged_modulus")),
  _a0(declareProperty<Real>("a0")),
  _a1(declareProperty<Real>("a1")),
  _a2(declareProperty<Real>("a2")),
  _a3(declareProperty<Real>("a3")),
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
  _fluid_viscosity(declareProperty<Real>("viscosity_fluid")),
  _C_g(declareProperty<Real>("C_g")),
  _m1(declareProperty<Real>("m1")),
  _m2(declareProperty<Real>("m2")), 
  _structural_stress_coefficient(declareProperty<Real>("structural_stress_coefficient")),
  _grad_alpha_damagedvar(declareProperty<RealGradient>("gradient_alpha_damagedvar")),
  _grad_alpha_damagedvar_xdir(declareProperty<Real>("gradient_alpha_damagedvar_xdir")),
  _grad_alpha_damagedvar_ydir(declareProperty<Real>("gradient_alpha_damagedvar_ydir")),
  _velgrad_L(declareProperty<RankTwoTensor>("velgrad_L")),
  //--------------------------------------------------------------//
  //input values
  _lambda_o_value(getParam<Real>("lambda_o")),
  _shear_modulus_o_value(getParam<Real>("shear_modulus_o")),
  _xi_0_value(getParam<Real>("xi_0")),
  _xi_d_value(getParam<Real>("xi_d")),
  _chi_value(getParam<Real>("chi")),
  _C_g_value(getParam<Real>("C_g")),
  _m1_value(getParam<Real>("m1")),
  _m2_value(getParam<Real>("m2")),
  _permeability_solid_o(getParam<Real>("permeability_solid_o")),
  _solid_bulk_modulus_s(getParam<Real>("solid_bulk_modulus_s")),
  _solid_bulk_modulus_g(getParam<Real>("solid_bulk_modulus_g")),
  _fluid_bulk_modulus(getParam<Real>("fluid_bulk_modulus")),
  _porosity_solid_o(getParam<Real>("porosity_solid_o")),
  _initial_viscosity_fluid(getParam<Real>("initial_viscosity_fluid")),
  _b(getParam<Real>("permeability_evolution_with_damage")),
  _DHo(getParam<Real>("initial_grain_size")),
  _DHu(getParam<Real>("ultimate_grain_size")),
  _alpha_damagedvar_aux(coupledValue("alpha_damagedvar_aux")),
  _B_damagedvar_aux(coupledValue("B_damagedvar_aux")),
  _structural_stress_coefficient_aux(coupledValue("structural_stress_coefficient")),
  _grad_alpha_damagedvar_value(coupledGradient("alpha_damagedvar_aux")),
  //--------------------------------------------------------------//
  _grad_vel_x(coupledGradient("vel_x")),
  _grad_vel_y(coupledGradient("vel_y")),
  _grad_vel_z(coupledGradient("vel_z")),
  //---------------------------------------------------------------//
  //Invariants, Jp and Dp - USE OLD VALUES TO BREAK CYCLIC DEPENDENCY
  _I1(getMaterialPropertyOldByName<Real>(getParam<std::string>("base_name") + "first_elastic_strain_invariant")),
  _Fe(getMaterialPropertyOldByName<RankTwoTensor>(getParam<std::string>("base_name") + "elastic_deformation_gradient")),
  _xi(getMaterialPropertyOldByName<Real>(getParam<std::string>("base_name") + "strain_invariant_ratio")),
  _Jp(getMaterialPropertyOldByName<Real>(getParam<std::string>("base_name") + "plastic_deformation_gradient_det")),
  _Dp(getMaterialPropertyOldByName<RankTwoTensor>(getParam<std::string>("base_name") + "plastic_strain_rate")),
  _Dp_old(getMaterialPropertyOlderByName<RankTwoTensor>(getParam<std::string>("base_name") + "plastic_strain_rate")),
  _phi_p_old(getMaterialPropertyOldByName<Real>(getParam<std::string>("base_name") + "plastic_porosity")),
  //---------------------------------------------------------------//
  //use spatial cg
  _use_spatial_cg(getParam<bool>("use_spatial_cg")),
  _cg_aux(_use_spatial_cg ? coupledValue("cg_aux") : _zero),
  //---------------------------------------------------------------//
  //use state dependent variables
  _use_state_var_evolution(getParam<bool>("use_state_var_evolution")),
  _const_A(getParam<Real>("const_A")),
  _const_B(getParam<Real>("const_B")),
  _const_theta_o(getParam<Real>("const_theta_o")),
  _initial_theta0(getParam<Real>("initial_theta0")),
  _use_state_var_evolution_mat(declareProperty<bool>("use_state_var_evolution_mat")),
  _const_A_mat(declareProperty<Real>("const_A_mat")),
  _const_B_mat(declareProperty<Real>("const_B_mat")),
  _const_theta_o_mat(declareProperty<Real>("const_theta_o_mat")),
  _initial_theta0_mat(declareProperty<Real>("initial_theta0_mat"))
  //---------------------------------------------------------------//
{
}

//Rules:See https://github.com/idaholab/moose/discussions/19450
//Only the object that declares the material property can assign values to it.
//Objects can request material properties, gaining read-only access to their values.
//When any object (including the object that declares it) requests the old value of a material property, that property becomes "stateful".
//All stateful material properties must be initialized within the initQpStatefulProperties call. 
//
void 
DiffusedDamageBreakageMaterialMainApp::initQpStatefulProperties()
{
  /* compute _gamma_damaged_r_mat */
  computegammar();

  /* update damage variable and breakage variable */
  updatedamagebreakage();

  /* compute modulus: _lambda, _shear_modulus, _damaged_modulus */
  updatemodulus();

  /* compute coefficients: a0 a1 a2 a3 */
  computecoefficients();

  /* compute: biot moduli and coefficients and permeabilities */
  updateporosolid();
  updateporogranular();

  /* compute L matrix */
  buildLmatrix();

  /* use state evolution */
  if (_use_state_var_evolution){ usestatevar();}

  /* compute spatial gradient of damage variable */
  _structural_stress_coefficient[_qp] = _structural_stress_coefficient_aux[_qp];
  _grad_alpha_damagedvar[_qp] = _grad_alpha_damagedvar_value[_qp];

  /* compute constant material properties: Cg, m1, m2 */
  if (_use_spatial_cg){ acceptspatialCg();}
  else{_C_g[_qp] = _C_g_value;}

  _m1[_qp] = _m1_value;
  _m2[_qp] = _m2_value;

  _phi_p[_qp] = 0.0;
}

void
DiffusedDamageBreakageMaterialMainApp::computeQpProperties()
{
  /* compute _gamma_damaged_r_mat */
  computegammar();

  /* update damage variable and breakage variable */
  updatedamagebreakage();

  /* compute modulus: _lambda, _shear_modulus, _damaged_modulus */
  updatemodulus();

  /* compute coefficients: a0 a1 a2 a3 */
  computecoefficients();

  /* compute: biot moduli and coefficients and permeabilities */
  updateporosolid();
  updateporogranular();

  /* compute L matrix */
  buildLmatrix();

  /* use state evolution */
  if (_use_state_var_evolution){ usestatevar();}

  /* compute spatial gradient of damage variable */
  _structural_stress_coefficient[_qp] = _structural_stress_coefficient_aux[_qp];
  _grad_alpha_damagedvar[_qp] = _grad_alpha_damagedvar_value[_qp];

  //for debugging
  _grad_alpha_damagedvar_xdir[_qp] = _grad_alpha_damagedvar_value[_qp](0);
  _grad_alpha_damagedvar_ydir[_qp] = _grad_alpha_damagedvar_value[_qp](1);

  /* compute constant material properties: Cg, m1, m2 */
  if (_use_spatial_cg){ acceptspatialCg();}
  else{_C_g[_qp] = _C_g_value;}

  _m1[_qp] = _m1_value;
  _m2[_qp] = _m2_value;
}

void 
DiffusedDamageBreakageMaterialMainApp::computegammar()
{

  // Calculate each part of the expression
  Real term1 = -_xi_0_value * (-_lambda_o_value * pow(_xi_0_value, 2) + 6 * _lambda_o_value + 2 * _shear_modulus_o_value);
  Real term2_sqrt = sqrt((_lambda_o_value * pow(_xi_0_value, 2) + 2 * _shear_modulus_o_value) * 
                            (_lambda_o_value * pow(_xi_0_value, 4) - 12 * _lambda_o_value * pow(_xi_0_value, 2) + 36 * _lambda_o_value
                            - 6 * _shear_modulus_o_value * pow(_xi_0_value, 2) + 24 * _shear_modulus_o_value));
  Real denominator = 2 * (pow(_xi_0_value, 2) - 3);
  
  // Calculate gamma_r
  Real gamma_r = (term1 - term2_sqrt) / denominator;
  
  //save
  _gamma_damaged_r[_qp] = gamma_r;
}

void
DiffusedDamageBreakageMaterialMainApp::updatedamagebreakage()
{
  _alpha_damagedvar[_qp] = _alpha_damagedvar_aux[_qp];
  _B_damagedvar[_qp] = _B_damagedvar_aux[_qp];
}

void 
DiffusedDamageBreakageMaterialMainApp::updatemodulus()
{
  Real shear_modulus = _shear_modulus_o_value +  _alpha_damagedvar[_qp] * _xi_0_value * _gamma_damaged_r[_qp];
  Real gamma_damaged =  _alpha_damagedvar[_qp] * _gamma_damaged_r[_qp];

  _lambda[_qp] = _lambda_o_value;
  _shear_modulus[_qp] = shear_modulus;
  _damaged_modulus[_qp] = gamma_damaged;
}

void
DiffusedDamageBreakageMaterialMainApp::computecoefficients()
{
  //compute xi_1
  Real _xi_1 = _xi_0_value + sqrt( pow(_xi_0_value , 2) + 2 * _shear_modulus_o_value / _lambda_o_value );

  //compute alpha_cr | xi = 0
  Real alpha_cr_xi0 = alphacr_root1(0);

  //compute mu_cr
  Real mu_cr = _shear_modulus_o_value + alpha_cr_xi0 * _xi_0_value * _gamma_damaged_r[_qp];

  //a0
  Real a0 = _chi_value * mu_cr;

  //a1
  Real numerator_a1 = -2 * _chi_value * mu_cr * pow(_xi_1, 3) + 6 * _chi_value * mu_cr * _xi_1 * pow(_xi_d_value, 2) - 4 * _chi_value * mu_cr * pow(_xi_d_value, 3)
                      - 2 * _gamma_damaged_r[_qp] * pow(_xi_1, 3) * _xi_d_value + 2 * _gamma_damaged_r[_qp] * pow(_xi_1, 3) * _xi_0_value
                      + _lambda_o_value * pow(_xi_1, 3) * pow(_xi_d_value, 2) + 2 * _shear_modulus_o_value * pow(_xi_1, 3);
  Real denominator_a1 = 2 * pow(_xi_1, 3) * _xi_d_value - 4 * pow(_xi_1, 2) * pow(_xi_d_value, 2) + 2 * _xi_1 * pow(_xi_d_value, 3);
  Real a1 = numerator_a1 / denominator_a1;

  //a2
  Real numerator_a2 = 2 * _chi_value * mu_cr * pow(_xi_1, 3) - 3 * _chi_value * mu_cr * pow(_xi_1, 2) * _xi_d_value + _chi_value * mu_cr * pow(_xi_d_value, 3)
                       + 2 * _gamma_damaged_r[_qp] * pow(_xi_1, 3) * _xi_d_value - 2 * _gamma_damaged_r[_qp] * pow(_xi_1, 3) * _xi_0_value
                       - _lambda_o_value * pow(_xi_1, 3) * pow(_xi_d_value, 2) - 2 * _shear_modulus_o_value * pow(_xi_1, 3);
  Real denominator_a2 = pow(_xi_1, 4) * _xi_d_value - 2 * pow(_xi_1, 3) * pow(_xi_d_value, 2) + pow(_xi_1, 2) * pow(_xi_d_value, 3); 
  Real a2 = numerator_a2 / denominator_a2; 

  //a3
  Real numerator_a3 = -2 * _chi_value * mu_cr * pow(_xi_1, 2) + 4 * _chi_value * mu_cr * _xi_1 * _xi_d_value - 2 * _chi_value * mu_cr * pow(_xi_d_value, 2)
                       - 2 * _gamma_damaged_r[_qp] * pow(_xi_1, 2) * _xi_d_value + 2 * _gamma_damaged_r[_qp] * pow(_xi_1, 2) * _xi_0_value
                       + _lambda_o_value * pow(_xi_1, 2) * pow(_xi_d_value, 2) + 2 * _shear_modulus_o_value * pow(_xi_1, 2);
  Real denominator_a3 = 2 * pow(_xi_1, 4) * _xi_d_value - 4 * pow(_xi_1, 3) * pow(_xi_d_value, 2) + 2 * pow(_xi_1, 2) * pow(_xi_d_value, 3);
  Real a3 = numerator_a3 / denominator_a3; 

  //save
  _a0[_qp] = a0;
  _a1[_qp] = a1;
  _a2[_qp] = a2;
  _a3[_qp] = a3;

}


void
DiffusedDamageBreakageMaterialMainApp::updateporosolid()
{
  // Get strain invariants
  Real I1 = _I1[_qp];
  Real xi = _xi[_qp];
  
  Real Je = _Fe[_qp].det();

  // Get damage/breakage variables
  Real alpha_damage = _alpha_damagedvar[_qp];
  Real B_breakage = _B_damagedvar[_qp];

  // Solid bulk modulus (constant for solid grains)
  Real K_s = _solid_bulk_modulus_s;
  
  // Fluid bulk modulus
  Real K_f = _fluid_bulk_modulus;
  
  // Compute drained bulk modulus K_d
  Real K_d = _lambda[_qp] + (2.0/3.0) * _shear_modulus[_qp] - (2.0/3.0) * _damaged_modulus[_qp] * xi;
  
  // Compute Biot coefficient for solid phase
  Real alpha_s = 1.0 - K_d/K_s;
  
  // Compute porosity evolution for solid phase
  Real porosity_s = 1 - (1 - _porosity_solid_o) * exp(-I1);
  // Real porosity_s = 1 - (1 - _porosity_solid_o) / Je;
  
   // Compute Biot modulus for solid phase
  Real one_over_Storage = (K_s*K_f)/(porosity_s * K_f + (alpha_s - porosity_s) * K_s);
  
  // Compute permeability for solid phase
  Real perm_s = _permeability_solid_o * pow(porosity_s/_porosity_solid_o, 3.0) * exp(_b * alpha_damage);
  
  // Save solid phase properties
  _Biot_coeff_s[_qp] = alpha_s;
  _Biot_modulus_s[_qp] = one_over_Storage;
  _perm_s[_qp] = (1 - B_breakage) * perm_s/_initial_viscosity_fluid;

  // Compute critical damage value
  Real alpha_cr = alphacr_root1(xi);

  // Determine critical porosity and permeability based on damage state
  const Real tolerance = 1e-10;  // Small tolerance for floating point comparison

  if (std::abs(_alpha_damagedvar[_qp] - alpha_cr) < tolerance) {
    // At critical damage: use current solid phase properties
    _phi_cr[_qp] = porosity_s;
    _perm_cr[_qp] = perm_s;
  } else {
    // Below critical damage: use initial solid properties
    _phi_cr[_qp] = _porosity_solid_o;
    _perm_cr[_qp] = _permeability_solid_o;
  }
  
}

void
DiffusedDamageBreakageMaterialMainApp::updateporogranular()
{
  // Get strain invariants 
  Real I1 = _I1[_qp];
  Real xi = _xi[_qp];

  Real Je = _Fe[_qp].det();

  // Get damage/breakage variables
  Real alpha_damage = _alpha_damagedvar[_qp];
  Real B_breakage = _B_damagedvar[_qp];

  // Fluid bulk modulus
  Real K_f = _fluid_bulk_modulus;

  // Solid bulk modulus (constant for solid grains)
  Real K_s_solid = _solid_bulk_modulus_s;

  // Solid bulk modulus after crushing
  Real K_s_crushed = _solid_bulk_modulus_g;
  
  // Compute bulk modulus for granular phase
  Real K_d_granular = 2 * _a2[_qp] + _a3[_qp] * (6.0 - (4.0/3.0) *  xi * xi ) * xi 
               + (2.0/3.0) * _a0[_qp] + (2.0/3.0) * _a1[_qp] * xi;

  // Solid bulk modulus for granular material 
  Real K_s_granular = (1 - B_breakage) * K_s_solid + B_breakage * K_s_crushed;

  // Compute Biot coefficient for granular phase
  Real alpha_g = 1.0 - K_d_granular/K_s_granular;

  // Compute elastic porosity evolution 
  Real porosity_e = 1 - (1 - _porosity_solid_o) * exp(-I1);
  // Real porosity_e = 1 - (1 - _porosity_solid_o) / Je;

  // Compute plastic porosity evolution 
  Real dporosity_pdt_new = (1 - _phi_p[_qp]) * _Dp[_qp].trace();
  Real dporosity_pdt_old = (1 - _phi_p_old[_qp]) * _Dp_old[_qp].trace();

  Real porosity_p = _phi_p_old[_qp] + 0.5 * _dt * (dporosity_pdt_old + dporosity_pdt_new);

  // Compute total porosity
  Real porosity = porosity_e + porosity_p;

  // Compute Biot modulus for solid phase
  Real one_over_Storage = (K_s_granular*K_f)/(porosity * K_f + (alpha_g - porosity) * K_s_granular);

  // Compute harmonic mean grain size for current distribution
  Real DH = (1 - B_breakage) * _DHo + B_breakage * _DHu;

 // Compute permeability for solid phase
  // k = k_cr * (φ/φ_0)^n * (DH/DH_0)^2  where n is typically 3
  // Real perm_g = _perm_cr[_qp] * pow(porosity/_phi_cr[_qp], 3.0) * pow(DH/_DHo, 2.0);
  Real perm_g = _perm_cr[_qp] * pow(porosity/_phi_cr[_qp], 3.0);

  // Save granular phase properties
  _Biot_coeff_g[_qp] = alpha_g;
  _Biot_modulus_g[_qp] = one_over_Storage;
  _phi_p[_qp] = porosity_p;
  _perm_g[_qp] = B_breakage * perm_g /_initial_viscosity_fluid;

  Real term22 = (1 - B_breakage) * _Biot_coeff_s[_qp] * _Biot_modulus_s[_qp] + B_breakage * _Biot_coeff_g[_qp] * _Biot_modulus_g[_qp];
  Real term33 = (1 - B_breakage) * _Biot_modulus_s[_qp] + B_breakage * _Biot_modulus_g[_qp];

  _Biot_modulus_eff[_qp] = term33;
  _fluid_solid_coupling[_qp] = term22 / term33;
  _biot_coeff_eff[_qp] = term22 / term33;
  _fluid_viscosity[_qp] = _initial_viscosity_fluid;

}


// Function for alpha_func_root1
Real 
DiffusedDamageBreakageMaterialMainApp::alphacr_root1(Real xi) {

  Real term1 = _lambda_o_value * pow(xi, 3) - 6 * _lambda_o_value * _xi_0_value + 6 * _shear_modulus_o_value * xi - 8 * _shear_modulus_o_value * _xi_0_value;
  Real term2 = std::sqrt(_lambda_o_value * _lambda_o_value * pow(xi, 6) 
                            - 12 * _lambda_o_value * _lambda_o_value * pow(xi, 3) * _xi_0_value 
                            + 36 * _lambda_o_value * _lambda_o_value * _xi_0_value * _xi_0_value 
                            + 12 * _lambda_o_value * _shear_modulus_o_value * pow(xi, 4) 
                            - 16 * _lambda_o_value * _shear_modulus_o_value * pow(xi, 3) * _xi_0_value 
                            - 72 * _lambda_o_value * _shear_modulus_o_value * pow(xi, 2) 
                            + 72 * _lambda_o_value * _shear_modulus_o_value * xi * _xi_0_value 
                            + 72 * _lambda_o_value * _shear_modulus_o_value 
                            - 12 * _shear_modulus_o_value * _shear_modulus_o_value * pow(xi, 2) 
                            + 48 * _shear_modulus_o_value * _shear_modulus_o_value);
  Real denominator = 2 * _gamma_damaged_r[_qp] * (3 * pow(xi, 2) - 6 * xi * _xi_0_value + 4 * _xi_0_value * _xi_0_value - 3);
  return (term1 - term2) / denominator;
}

//build the L matrix
//L matrix is the gradient of velocity
void
DiffusedDamageBreakageMaterialMainApp::buildLmatrix()
{
  _velgrad_L[_qp](0,0) = (_grad_vel_x)[_qp](0); _velgrad_L[_qp](0,1) = (_grad_vel_x)[_qp](1); _velgrad_L[_qp](0,2) = (_grad_vel_x)[_qp](2);
  _velgrad_L[_qp](1,0) = (_grad_vel_y)[_qp](0); _velgrad_L[_qp](1,1) = (_grad_vel_y)[_qp](1); _velgrad_L[_qp](1,2) = (_grad_vel_y)[_qp](2);
  _velgrad_L[_qp](2,0) = (_grad_vel_z)[_qp](0); _velgrad_L[_qp](2,1) = (_grad_vel_z)[_qp](1); _velgrad_L[_qp](2,2) = (_grad_vel_z)[_qp](2);
}

void 
DiffusedDamageBreakageMaterialMainApp::acceptspatialCg()
{
  _C_g[_qp] = _cg_aux[_qp];
}

void 
DiffusedDamageBreakageMaterialMainApp::usestatevar()
{
  _use_state_var_evolution_mat[_qp] = _use_state_var_evolution;
  _const_A_mat[_qp] = _const_A;
  _const_B_mat[_qp] = _const_B;
  _const_theta_o_mat[_qp] = _const_theta_o;
  _initial_theta0_mat[_qp] = _initial_theta0;
}