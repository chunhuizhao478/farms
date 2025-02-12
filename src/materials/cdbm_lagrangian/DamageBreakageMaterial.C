//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DamageBreakageMaterial.h"

/**
 *  Material used in damage-breakage large deformation formulation
 *  Created by Chunhui Zhao, Aug 5th, 2024
 */
registerMooseObject("farmsApp", DamageBreakageMaterial);

InputParameters
DamageBreakageMaterial::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Material used in three field poro dynamics simulations");
  params.addRequiredParam<Real>(        "lambda_o", "initial lambda constant value");
  params.addRequiredParam<Real>( "shear_modulus_o", "initial shear modulus value");
  params.addRequiredParam<Real>(            "xi_0", "strain invariants ratio: onset of damage evolution");
  params.addRequiredParam<Real>(            "xi_d", "strain invariants ratio: onset of breakage healing");
  params.addRequiredParam<Real>(          "xi_min", "strain invariants ratio: minimum allowable value");
  params.addRequiredParam<Real>(          "xi_max", "strain invariants ratio: maximum allowable value");
  params.addRequiredParam<Real>(             "chi", "coefficient of energy ratio Fb/Fs = chi < 1");
  params.addParam<Real>(     "Cd_constant", -1, "coefficient gives positive damage evolution");
  params.addParam<Real>(             "C_1", -1, "coefficient of healing for damage evolution");
  params.addParam<Real>(             "C_2", -1, "coefficient of healing for damage evolution");
  params.addRequiredParam<Real>(      "beta_width", "coefficient gives width of transitional region");
  params.addParam<Real>( "CdCb_multiplier", -1, "multiplier between Cd and Cb");
  params.addParam<Real>(    "CBH_constant", -1, "constant CBH value");
  params.addParam<Real>(             "C_g", -1, "compliance or fluidity of the fine grain granular material");
  params.addRequiredParam<Real>(              "m1", "coefficient of power law indexes");
  params.addRequiredParam<Real>(              "m2", "coefficient of power law indexes");
  params.addParam<bool>("use_xi0_aux", false, "Whether to use xi_0 from aux variable");
  params.addCoupledVar("xi0_aux", "AuxVariable for xi_0");
  params.addParam<bool>("use_shear_modulus_o_aux", false, "Whether to use shear_modulus_o from aux variable");
  params.addCoupledVar("shear_modulus_o_aux", "AuxVariable for shear_modulus_o");
  params.addParam<bool>("use_nonlocal_xi", false, "Whether to use nonlocal xi variable");
  params.addCoupledVar("nonlocal_xi", "Nonlocal xi variable");
  params.addParam<bool>("use_cd_aux", false, "Whether to use Cd from aux variable");
  params.addCoupledVar("Cd_constant_aux", "AuxVariable for Cd block-restricted");
  params.addParam<bool>("use_cb_multiplier_aux", false, "Whether to use Cb multiplier from aux variable");
  params.addCoupledVar("Cb_multiplier_aux", "AuxVariable for Cb multiplier block-restricted");
  params.addParam<bool>("use_cbh_aux", false, "Whether to use CBH from aux variable");
  params.addCoupledVar("CBH_aux", "AuxVariable for CBH block-restricted");
  params.addParam<bool>("use_c1_aux", false, "Whether to use C1 from aux variable");
  params.addCoupledVar("C1_aux", "AuxVariable for C1 block-restricted");
  params.addParam<bool>("use_c2_aux", false, "Whether to use C2 from aux variable");
  params.addCoupledVar("C2_aux", "AuxVariable for C2 block-restricted");
  params.addParam<bool>("use_cg_aux", false, "Whether to use Cg from aux variable");
  params.addCoupledVar("Cg_aux", "AuxVariable for Cg block-restricted");
  params.addParam<bool>("use_const_xi_aux", false, "Whether to use xi_0 from aux variable");
  params.addCoupledVar("const_xi_aux", "AuxVariable for xi");
  params.addParam<int>("const_xi_block_id", -1, "Block ID to apply the constant xi, currently it has to be one block");
  params.addParam<bool>("use_cd_strain_dependent", false, "Whether to use strain-dependent Cd");
  params.addParam<Real>("m_exponent", -1, "Exponent for strain-dependent Cd");
  params.addParam<Real>("strain_rate_hat", -1, "Strain rate for strain-dependent Cd");
  params.addParam<Real>("cd_hat", -1, "Cd value for strain-dependent Cd");
  params.addParam<int>("straindep_block_id_applied", -1, "Block ID to apply the rate-dependent Cd, currently it has to be one block");
  params.addParam<bool>("use_total_strain_rate", false, "Whether to use total strain rate, default is to use elastic strain rate");
  // Add option to use pore pressure to decrease mean stress
  params.addParam<bool>("use_pore_pressure", false,
                        "Flag to use pore pressure to decrease mean stress");
  // Optional coupled variable: the user should supply the name of the pore pressure variable if used.
  params.addCoupledVar("pore_pressure", "Name of the pore pressure coupled variable");
  params.addParam<bool>("use_overstress", false,
                        "Flag to use overstress to nucleate the rupture");
  params.addCoupledVar("overstress", "Name of the pore pressure coupled variable");
  // Build initial damage profile inside this material object
  params.addParam<bool>("build_param_use_initial_damage_time_dependent_mat", false,
                        "Flag to build initial damage profile inside this material object");
  params.addParam<Real>("build_param_peak_value", -1,
                        "Peak value of the initial damage profile");
  params.addParam<Real>("build_param_sigma", -1,
                        "Sigma value of the initial damage profile");
  params.addParam<Real>("build_param_len_of_fault", -1,
                        "Length of the fault for the initial damage profile");
  // Build time dependent damage perturbation inside this material object
  params.addParam<bool>("perturbation_build_param_use_damage_perturb", false,
                        "Use damage perturbation for the initial damage profile");
  params.addParam<std::vector<Real>>("perturbation_build_param_nucl_center", std::vector<Real>(),
                        "Center of the nucleation for the damage perturbation");
  params.addParam<Real>("perturbation_build_param_length", -1,
                        "Length of the nucleation for the damage perturbation");
  params.addParam<Real>("perturbation_build_param_peak_value", -1,
                        "Peak value of the nucleation for the damage perturbation");
  params.addParam<Real>("perturbation_build_param_sigma", -1,
                        "Sigma value of the nucleation for the damage perturbation");
  params.addParam<Real>("perturbation_build_param_thickness", -1,
                        "Thickness of the nucleation for the damage perturbation");
  params.addParam<Real>("perturbation_build_param_duration", -1,
                        "Duration of the nucleation for the damage perturbation");  
  return params;
}

DamageBreakageMaterial::DamageBreakageMaterial(const InputParameters & parameters)
  : Material(parameters),
  _alpha_damagedvar(declareProperty<Real>("alpha_damagedvar")),
  _B_breakagevar(declareProperty<Real>("B_damagedvar")),
  _initial_damage_time_dependent_mat(declareProperty<Real>("initial_damage_time_dependent_material")),
  _damage_perturbation(declareProperty<Real>("damage_perturb_material")),
  _damage_perturbation_old(getMaterialPropertyOldByName<Real>("damage_perturb_material")),
  _lambda(declareProperty<Real>("lambda_const")),
  _shear_modulus(declareProperty<Real>("shear_modulus")),
  _damaged_modulus(declareProperty<Real>("damaged_modulus")),
  _gamma_damaged_r(declareProperty<Real>("gamma_damaged_r")),
  _a0(declareProperty<Real>("a0")),
  _a1(declareProperty<Real>("a1")),
  _a2(declareProperty<Real>("a2")),
  _a3(declareProperty<Real>("a3")),
  _C_g(declareProperty<Real>("C_g")),
  _m1(declareProperty<Real>("m1")),
  _m2(declareProperty<Real>("m2")),  
  _xi_1_mat(declareProperty<Real>("xi_1")),  
  _xi_0_mat(declareProperty<Real>("xi_0")),
  _shear_modulus_o_mat(declareProperty<Real>("shear_modulus_o_mat")),
  _Cd_rate_dependent(declareProperty<Real>("Cd_rate_dependent")),
  _strain_dir0_positive(declareProperty<Real>("strain_dir0_positive")),
  _strain_dir0_positive_old(getMaterialPropertyOldByName<Real>("strain_dir0_positive")),
  _alpha_damagedvar_old(getMaterialPropertyOldByName<Real>("alpha_damagedvar")),
  _B_breakagevar_old(getMaterialPropertyOldByName<Real>("B_damagedvar")),
  _I2_old(getMaterialPropertyOldByName<Real>("second_elastic_strain_invariant")),
  _xi_old(getMaterialPropertyOldByName<Real>("strain_invariant_ratio")),
  _initial_damage(getMaterialProperty<Real>("initial_damage")),
  _gamma_damaged_r_old(getMaterialPropertyOldByName<Real>("gamma_damaged_r")),
  _a0_old(getMaterialPropertyOldByName<Real>("a0")),
  _a1_old(getMaterialPropertyOldByName<Real>("a1")),
  _a2_old(getMaterialPropertyOldByName<Real>("a2")),
  _a3_old(getMaterialPropertyOldByName<Real>("a3")),
  _elastic_strain_old(getMaterialPropertyOldByName<RankTwoTensor>("green_lagrange_elastic_strain")),
  _total_lagrange_strain_old(getMaterialPropertyOldByName<RankTwoTensor>("total_lagrange_strain")),
  _lambda_o(getParam<Real>("lambda_o")),
  _shear_modulus_o(getParam<Real>("shear_modulus_o")),
  _xi_d(getParam<Real>("xi_d")),
  _xi_min(getParam<Real>("xi_min")),
  _xi_max(getParam<Real>("xi_max")),
  _chi(getParam<Real>("chi")),
  _Cd_constant(getParam<Real>("Cd_constant")), 
  _C1(getParam<Real>("C_1")),
  _C2(getParam<Real>("C_2")),
  _beta_width(getParam<Real>("beta_width")),
  _CdCb_multiplier(getParam<Real>("CdCb_multiplier")),
  _CBH_constant(getParam<Real>("CBH_constant")),
  _C_g_value(getParam<Real>("C_g")),
  _m1_value(getParam<Real>("m1")),
  _m2_value(getParam<Real>("m2")),
  _use_xi0_aux(getParam<bool>("use_xi0_aux")),
  _xi0_value(getParam<Real>("xi_0")),
  _xi0_aux(_use_xi0_aux ? &coupledValue("xi0_aux") : nullptr),
  _use_shear_modulus_o_aux(getParam<bool>("use_shear_modulus_o_aux")),
  _shear_modulus_o_value(getParam<Real>("shear_modulus_o")),
  _shear_modulus_o_aux(_use_shear_modulus_o_aux ? &coupledValue("shear_modulus_o_aux") : nullptr),
  _use_nonlocal_xi(getParam<bool>("use_nonlocal_xi")),
  _nonlocal_xi(_use_nonlocal_xi ? &coupledValueOld("nonlocal_xi") : nullptr),
  _use_cd_aux(getParam<bool>("use_cd_aux")),
  _cd_aux(_use_cd_aux ? &coupledValue("Cd_constant_aux") : nullptr),
  _use_cb_multiplier_aux(getParam<bool>("use_cb_multiplier_aux")),
  _cb_multiplier_aux(_use_cb_multiplier_aux ? &coupledValue("Cb_multiplier_aux") : nullptr),
  _use_cbh_aux(getParam<bool>("use_cbh_aux")),
  _cbh_aux(_use_cbh_aux ? &coupledValue("CBH_aux") : nullptr),
  _use_c1_aux(getParam<bool>("use_c1_aux")),
  _c1_aux(_use_c1_aux ? &coupledValue("C1_aux") : nullptr),
  _use_c2_aux(getParam<bool>("use_c2_aux")),
  _c2_aux(_use_c2_aux ? &coupledValue("C2_aux") : nullptr),
  _use_cg_aux(getParam<bool>("use_cg_aux")),
  _cg_aux(_use_cg_aux ? &coupledValue("Cg_aux") : nullptr),
  _use_const_xi_aux(getParam<bool>("use_const_xi_aux")),
  _const_xi_aux(_use_const_xi_aux ? &coupledValue("const_xi_aux") : nullptr),
  _const_xi_block_id(getParam<int>("const_xi_block_id")),
  _use_cd_strain_dependent(getParam<bool>("use_cd_strain_dependent")),
  _block_id(0),  // Initialize block ID
  _m_exponent(getParam<Real>("m_exponent")),
  _strain_rate_hat(getParam<Real>("strain_rate_hat")),
  _cd_hat(getParam<Real>("cd_hat")),
  _straindep_block_id_applied(getParam<int>("straindep_block_id_applied")),
  _use_total_strain_rate(getParam<bool>("use_total_strain_rate")),
  _use_pore_pressure(getParam<bool>("use_pore_pressure")),
  _pore_pressure(_use_pore_pressure ? &coupledValue("pore_pressure") : nullptr),
  _pore_pressure_mat(declareProperty<Real>("pore_pressure_mat")),
  _use_overstress(getParam<bool>("use_overstress")),
  _overstress(_use_overstress ? &coupledValue("overstress") : nullptr),
  _use_overstress_mat(declareProperty<bool>("use_overstress_mat")),
  _overstress_mat(declareProperty<Real>("overstress_mat")),
  // Build initial damage profile inside this material object
  _build_param_use_initial_damage_time_dependent_mat(getParam<bool>("build_param_use_initial_damage_time_dependent_mat")),
  _build_param_peak_value(getParam<Real>("build_param_peak_value")),
  _build_param_sigma(getParam<Real>("build_param_sigma")),
  _build_param_len_of_fault(getParam<Real>("build_param_len_of_fault")),
  // Build time dependent damage perturbation inside this material object
  _perturbation_build_param_use_damage_perturb(getParam<bool>("perturbation_build_param_use_damage_perturb")),
  _perturbation_build_param_nucl_center(getParam<std::vector<Real>>("perturbation_build_param_nucl_center")),
  _perturbation_build_param_length(getParam<Real>("perturbation_build_param_length")),
  _perturbation_build_param_peak_value(getParam<Real>("perturbation_build_param_peak_value")),
  _perturbation_build_param_sigma(getParam<Real>("perturbation_build_param_sigma")),
  _perturbation_build_param_thickness(getParam<Real>("perturbation_build_param_thickness")),
  _perturbation_build_param_duration(getParam<Real>("perturbation_build_param_duration"))
{
  if (_use_xi0_aux && !parameters.isParamSetByUser("xi0_aux"))
    mooseError("Must specify xi0_aux when use_xi0_aux = true");
  if (_use_shear_modulus_o_aux && !parameters.isParamSetByUser("shear_modulus_o_aux"))
    mooseError("Must specify shear_modulus_o_aux when use_shear_modulus_o_aux = true");
  if (_use_nonlocal_xi && !parameters.isParamSetByUser("nonlocal_xi"))
    mooseError("Must specify nonlocal_xi when use_nonlocal_xi = true");
  if (_use_cd_aux && !parameters.isParamSetByUser("Cd_constant_aux"))
    mooseError("Must specify Cd_constant_aux when use_cd_aux = true");
  if (_use_cb_multiplier_aux && !parameters.isParamSetByUser("Cb_multiplier_aux"))
    mooseError("Must specify Cb_multiplier_aux when use_cb_multiplier_aux = true");
  if (_use_cbh_aux && !parameters.isParamSetByUser("CBH_aux"))
    mooseError("Must specify CBH_aux when use_cbh_aux = true");
  if (_use_c1_aux && !parameters.isParamSetByUser("C1_aux"))
    mooseError("Must specify C1_aux when use_c1_aux = true");
  if (_use_c2_aux && !parameters.isParamSetByUser("C2_aux"))
    mooseError("Must specify C2_aux when use_c2_aux = true");
  if (_Cd_constant < 0 && !_use_cd_aux)
    mooseError("Cd_constant must be set to a positive value or use_cd_aux must be set to true");
  if (_C1 < 0 && !_use_c1_aux)
    mooseError("C1 must be set to a positive value or use_c1_aux must be set to true");
  if (_C2 < 0 && !_use_c2_aux)
    mooseError("C2 must be set to a positive value or use_c2_aux must be set to true");
  if (_CdCb_multiplier < 0 && !_use_cb_multiplier_aux)
    mooseError("CdCb_multiplier must be set to a positive value or use_cb_multiplier_aux must be set to true");
  if (_CBH_constant < 0 && !_use_cbh_aux)
    mooseError("CBH_constant must be set to a positive value or use_cbh_aux must be set to true");
  if (_Cd_constant > 0 && _use_cd_aux)
    mooseError("Global Param Cd_constant must not be set when use_cd_aux is set to true");
  if (_C1 > 0 && _use_c1_aux)
    mooseError("Global Param C1 must not be set when use_c1_aux is set to true");
  if (_C2 > 0 && _use_c2_aux)
    mooseError("Global Param C2 must not be set when use_c2_aux is set to true");
  if (_CdCb_multiplier > 0 && _use_cb_multiplier_aux)
    mooseError("Global Param CdCb_multiplier must not be set when use_cb_multiplier_aux is set to true");
  if (_CBH_constant > 0 && _use_cbh_aux)
    mooseError("Global Param CBH_constant must not be set when use_cbh_aux is set to true");
  if ((_use_cd_strain_dependent) && (_m_exponent < 0))
    mooseError("m_exponent must be set to a positive value when use_cd_strain_dependent is set to true");
  if ((_use_cd_strain_dependent) && (_strain_rate_hat < 0))
    mooseError("strain_rate_hat must be set to a positive value when use_cd_strain_dependent is set to true");
  if ((_use_cd_strain_dependent) && (_cd_hat < 0))
    mooseError("cd_hat must be set to a positive value when use_cd_strain_dependent is set to true");
  if ((_use_cd_strain_dependent) && (_straindep_block_id_applied < 0))
    mooseError("block_id_applied must be set to a positive value when use_cd_strain_dependent is set to true");
  if ((_use_const_xi_aux) && (_const_xi_block_id < 0))
    mooseError("const_xi_block_id must be set to a positive value when use_const_xi_aux is set to true");
  // If use_pore_pressure is true but no pore_pressure is provided, you may throw an error:
  if (_use_pore_pressure && !parameters.isParamSetByUser("pore_pressure"))
    mooseError("Must specify pore_pressure when use_pore_pressure = true");  
  if (_use_overstress && !parameters.isParamSetByUser("overstress"))
    mooseError("Must specify overstress when use_overstress = true");
  // Check parameters for building initial damage
  if (_build_param_use_initial_damage_time_dependent_mat && (_build_param_peak_value < 0))
    mooseError("build_param_peak_value must be set to a positive value when build_param_use_initial_damage_time_dependent_mat is set to true");
  if (_build_param_use_initial_damage_time_dependent_mat && (_build_param_sigma < 0))
    mooseError("build_param_sigma must be set to a positive value when build_param_use_initial_damage_time_dependent_mat is set to true");
  if (_build_param_use_initial_damage_time_dependent_mat && (_build_param_len_of_fault < 0))
    mooseError("build_param_len_of_fault must be set to a positive value when build_param_use_initial_damage_time_dependent_mat is set to true");
  // Check parameters for building time dependent damage perturbation
  if (_perturbation_build_param_use_damage_perturb && (_perturbation_build_param_nucl_center.size() != 2))
    mooseError("perturbation_build_param_nucl_center must be set to a 2D vector when perturbation_build_param_use_damage_perturb is set to true");
  if (_perturbation_build_param_use_damage_perturb && (_perturbation_build_param_length < 0))
    mooseError("perturbation_build_param_length must be set to a positive value when perturbation_build_param_use_damage_perturb is set to true");
  if (_perturbation_build_param_use_damage_perturb && (_perturbation_build_param_peak_value < 0))
    mooseError("perturbation_build_param_peak_value must be set to a positive value when perturbation_build_param_use_damage_perturb is set to true");    
  if (_perturbation_build_param_use_damage_perturb && (_perturbation_build_param_sigma < 0))
    mooseError("perturbation_build_param_sigma must be set to a positive value when perturbation_build_param_use_damage_perturb is set to true");
  if (_perturbation_build_param_use_damage_perturb && (_perturbation_build_param_thickness < 0))
    mooseError("perturbation_build_param_thickness must be set to a positive value when perturbation_build_param_use_damage_perturb is set to true");
  if (_perturbation_build_param_use_damage_perturb && (_perturbation_build_param_duration < 0))
    mooseError("perturbation_build_param_duration must be set to a positive value when perturbation_build_param_use_damage_perturb is set to true");
}

//Rules:See https://github.com/idaholab/moose/discussions/19450
//Only the object that declares the material property can assign values to it.
//Objects can request material properties, gaining read-only access to their values.
//When any object (including the object that declares it) requests the old value of a material property, that property becomes "stateful".
//All stateful material properties must be initialized within the initQpStatefulProperties call. 
//
void 
DamageBreakageMaterial::initQpStatefulProperties()
{

  // Get xi_0 value
  const Real _xi_0 = _use_xi0_aux ? (*_xi0_aux)[_qp] : _xi0_value;
  // Get shear_modulus_o value
  const Real _shear_modulus_o = _use_shear_modulus_o_aux ? (*_shear_modulus_o_aux)[_qp] : _shear_modulus_o_value;

  /* compute gamma_r */
  computegammar();

  /* compute a0 a1 a2 a3 coefficients */
  computecoefficients();

  _alpha_damagedvar[_qp] = _initial_damage[_qp]; //
  _B_breakagevar[_qp] = 0.0;

  _lambda[_qp] = _lambda_o;
  _shear_modulus[_qp] = _shear_modulus_o + _initial_damage[_qp] * _xi_0 * _gamma_damaged_r[_qp];
  _damaged_modulus[_qp] = _initial_damage[_qp] * _gamma_damaged_r[_qp];
  
  /* define Cg m1 m2 */
  _m1[_qp] = _m1_value;
  _m2[_qp] = _m2_value;

  // Get cg value
  _C_g[_qp] = _C_g_value;

  //initialize maximum principal strain 
  _strain_dir0_positive[_qp] = 0.0;

  //initialize pore pressure
  _pore_pressure_mat[_qp] = 0.0;

  //initialize overstress
  _use_overstress_mat[_qp] = false;
  _overstress_mat[_qp] = 0.0;

  //initialize initial damage time dependent material
  _initial_damage_time_dependent_mat[_qp] = 0.0;

  //initialize damage perturbation
  _damage_perturbation[_qp] = 0.0;

}

void
DamageBreakageMaterial::computeQpProperties()
{

  /* compute gamma_r */
  computegammar();

  /* compute a0 a1 a2 a3 coefficients */
  computecoefficients();

  // Build time dependent damage perturbation inside this material object
  if (_perturbation_build_param_use_damage_perturb){
    computedamageperturbation2D();
  }

  // Build initial damage profile inside this material object
  if (_build_param_use_initial_damage_time_dependent_mat){
    computeinitialdamage2D();
  }

  /* compute alpha at t_{n+1} using quantities from t_{n} */
  Real alpha_updated = updatedamage();

  /* compute B at t_{n+1} using quantities from t_{n} */
  updatebreakage();

  /* compute modulus at t_{n+1} using alpha at t_{n} */
  updatemodulus(alpha_updated);

  /* define Cg m1 m2 */
  _m1[_qp] = _m1_value;
  _m2[_qp] = _m2_value;

  // Get cg value
  _C_g[_qp] = _C_g_value;

  // Get xi_0 value
  const Real _xi_0 = _use_xi0_aux ? (*_xi0_aux)[_qp] : _xi0_value;
  _xi_0_mat[_qp] = _xi_0;
  // Get shear_modulus_o value
  const Real _shear_modulus_o = _use_shear_modulus_o_aux ? (*_shear_modulus_o_aux)[_qp] : _shear_modulus_o_value;
  _shear_modulus_o_mat[_qp] = _shear_modulus_o;

  // Get pore pressure value
  if (_use_pore_pressure)
    _pore_pressure_mat[_qp] = _use_pore_pressure && _pore_pressure ? (*_pore_pressure)[_qp] : 0.0;

  // Get overstress value
  if (_use_overstress)
    _use_overstress_mat[_qp] = true;
    _overstress_mat[_qp] = _use_overstress && _overstress ? (*_overstress)[_qp] : 0.0;

}

Real 
DamageBreakageMaterial::updatedamage()
{

  // Get the block ID for the current element
  _block_id = _current_elem->subdomain_id();

  // Get xi_0 value
  const Real _xi_0 = _use_xi0_aux ? (*_xi0_aux)[_qp] : _xi0_value;

  // Get xi value based on option
  Real xi = _use_nonlocal_xi ? (*_nonlocal_xi)[_qp] : _xi_old[_qp];

  //---------------------------buffer zone-------------------------------//
  // Get xi value based on whether using aux variable
  // this is block-restricted operation, only apply to buffer-zone block
  if ((_use_const_xi_aux) && (_block_id == _const_xi_block_id)){
    xi = (*_const_xi_aux)[_qp];
  }
  //---------------------------buffer zone-------------------------------//

  // Get Cd value based on option
  Real Cd = _use_cd_aux ? (*_cd_aux)[_qp] : _Cd_constant;

  // Get Cd value if using strain dependent
  RealVectorValue strain_in_crack_dir;
  const Real strain_dir0_positive_old = _strain_dir0_positive_old[_qp];

  // std::cout << "Block ID: " << _block_id << std::endl;
  if ((_use_cd_strain_dependent) && (_block_id == _straindep_block_id_applied )){ //avoid modifying outer block zero Cd
    //compute principal strain
    computePrincipalStrainAndOrientation(strain_in_crack_dir);
    //compute maximum tensile strain
    Real strain_dir0_positive = std::max(strain_in_crack_dir(0), 0.0); //maximum tensile strain
    //compute strain rate
    Real strain_rate = (strain_dir0_positive - strain_dir0_positive_old) / _dt; //strain rate
    //save strain_dir0_positive
    _strain_dir0_positive[_qp] = strain_dir0_positive;
    //compute rate-dependent Cd
    if (strain_rate < _strain_rate_hat){ strain_rate = _strain_rate_hat; } //avoid zero strain rate
    Cd = pow(10, 1 + _m_exponent * log(strain_rate/_strain_rate_hat)) * _cd_hat; //scale the Cd value //take constant Cd as minimum
  } //note: {} is needed for statements, otherwise weird error will occur

  // Save rate-denpendent Cd
  _Cd_rate_dependent[_qp] = Cd;

  // Get C1 value based on option
  const Real C1 = _use_c1_aux ? (*_c1_aux)[_qp] : _C1;

  // Get C2 value based on option
  Real C2 = _use_c2_aux ? (*_c2_aux)[_qp] : _C2;

  // Set initial C2 as nonzero value (for the use of aux variable)
  if ( C2 == 0 ){ C2 = 0.05; } 

  //compute forcing term
  Real alpha_forcingterm;
  if ( xi >= _xi_0 ){
    alpha_forcingterm = (1 - _B_breakagevar_old[_qp]) * ( Cd * _I2_old[_qp] * ( xi - _xi_0 ) );
  }
  else if ( xi < _xi_0 ){
    alpha_forcingterm = (1 - _B_breakagevar_old[_qp]) * ( C1 * std::exp(_alpha_damagedvar_old[_qp]/C2) * _I2_old[_qp] * ( xi - _xi_0 ) );
  }
  // else{
  //   mooseError("xi_old is OUT-OF-RANGE!.");   
  // }

  //update alpha at current time
  Real alpha_damagedvar = _alpha_damagedvar_old[_qp] + _dt * alpha_forcingterm;

  //check alpha within range
  if ( alpha_damagedvar < 0 ){ alpha_damagedvar = 0.0; }
  else if ( alpha_damagedvar > 1 ){ alpha_damagedvar = 1.0; }
  else{} 

  
  //Get currernt initial damage
  //This function will separate the case if time dependent initial damage is used
  Real get_correct_initial_damage = _initial_damage[_qp];
  if ( _build_param_use_initial_damage_time_dependent_mat ){
    get_correct_initial_damage = _initial_damage_time_dependent_mat[_qp]; 
  }

  //check below initial damage (fix initial damage)
  if ( alpha_damagedvar < get_correct_initial_damage ){ alpha_damagedvar = get_correct_initial_damage; }
  else{}

  _alpha_damagedvar[_qp] = alpha_damagedvar;

  return alpha_damagedvar;
}

void
DamageBreakageMaterial::updatebreakage()
{

  // Get xi_0 value
  const Real _xi_0 = _use_xi0_aux ? (*_xi0_aux)[_qp] : _xi0_value;
  // Get shear_modulus_o value
  const Real _shear_modulus_o = _use_shear_modulus_o_aux ? (*_shear_modulus_o_aux)[_qp] : _shear_modulus_o_value;

  // Get Cd value based on option
  Real Cd = _use_cd_aux ? (*_cd_aux)[_qp] : _Cd_constant;  
  // Get Cd value based on option
  const Real CdCb_multiplier = _use_cb_multiplier_aux ? (*_cb_multiplier_aux)[_qp] : _CdCb_multiplier;

  // Get Cbh value based on option
  const Real CBH_constant = _use_cbh_aux ? (*_cbh_aux)[_qp] : _CBH_constant;

  // Get Cd value if using strain dependent
  if (_use_cd_strain_dependent){
    Cd = _Cd_rate_dependent[_qp];
  }

  /* compute C_B based on C_d */
  Real C_B = CdCb_multiplier * Cd;

  //compute xi_1
  Real _xi_1 = _xi_0 + sqrt( pow(_xi_0 , 2) + 2 * _shear_modulus_o / _lambda_o );

  //save
  _xi_1_mat[_qp] = _xi_1;

  // Get xi value based on option
  Real xi = _use_nonlocal_xi ? (*_nonlocal_xi)[_qp] : _xi_old[_qp];  

  //----------------------------buffer zone-----------------------------//
  // Get the block ID for the current element
  _block_id = _current_elem->subdomain_id();
  // Get xi value based on whether using aux variable
  // this is block-restricted operation, only apply to buffer-zone block
  if ((_use_const_xi_aux) && (_block_id == _const_xi_block_id)){
    xi = (*_const_xi_aux)[_qp];  
  }
  //----------------------------buffer zone-----------------------------//

  //alphacr function
  Real alphacr;
  if ( xi < _xi_0 ){ alphacr = 1.0;} 
  else if ( xi > _xi_0 && xi <= _xi_1 ){ alphacr = alphacr_root1(xi);}
  else if ( xi > _xi_1 && xi <= _xi_max ){ alphacr = alphacr_root2(xi);}
  else{std::cout<<"xi: "<<xi<<std::endl; mooseError("xi exceeds the maximum allowable range!");}

  //compute forcing func
  Real Prob = 1.0 / ( std::exp( (alphacr - _alpha_damagedvar_old[_qp]) / _beta_width ) + 1.0 );
  Real B_forcingterm;
  if ( xi >= _xi_d ){
    B_forcingterm = 1.0 * C_B * Prob * (1-_B_breakagevar_old[_qp]) * _I2_old[_qp] * (xi - _xi_d); //could heal if xi < xi_0
  }
  else if ( xi < _xi_d ){
    B_forcingterm = 1.0 * CBH_constant * _I2_old[_qp] * ( xi - _xi_d );
  }
  // else{
  //   mooseError("xi_old is OUT-OF-RANGE!.");
  // }

  Real B_damagedvar = _B_breakagevar_old[_qp] + _dt * B_forcingterm;

  //check breakage within range
  if ( B_damagedvar < 0 ){ B_damagedvar = 0.0; }
  else if ( B_damagedvar > 1 ){ B_damagedvar = 1.0; }
  else{}   

  _B_breakagevar[_qp] = B_damagedvar;

}

void 
DamageBreakageMaterial::updatemodulus(Real alpha_updated)
{

  // Get xi_0 value
  const Real _xi_0 = _use_xi0_aux ? (*_xi0_aux)[_qp] : _xi0_value;
  // Get shear_modulus_o value
  const Real _shear_modulus_o = _use_shear_modulus_o_aux ? (*_shear_modulus_o_aux)[_qp] : _shear_modulus_o_value;

  Real shear_modulus = _shear_modulus_o + alpha_updated * _xi_0 * _gamma_damaged_r[_qp];
  Real gamma_damaged = alpha_updated * _gamma_damaged_r[_qp];

  _lambda[_qp] = _lambda_o;
  _shear_modulus[_qp] = shear_modulus;
  _damaged_modulus[_qp] = gamma_damaged;
}

void 
DamageBreakageMaterial::computegammar()
{

  // Get xi_0 value
  const Real _xi_0 = _use_xi0_aux ? (*_xi0_aux)[_qp] : _xi0_value;
  // Get shear_modulus_o value
  const Real _shear_modulus_o = _use_shear_modulus_o_aux ? (*_shear_modulus_o_aux)[_qp] : _shear_modulus_o_value;

  // Calculate each part of the expression
  Real term1 = -_xi_0 * (-_lambda_o * pow(_xi_0, 2) + 6 * _lambda_o + 2 * _shear_modulus_o);
  Real term2_sqrt = sqrt((_lambda_o * pow(_xi_0, 2) + 2 * _shear_modulus_o) * 
                            (_lambda_o * pow(_xi_0, 4) - 12 * _lambda_o * pow(_xi_0, 2) + 36 * _lambda_o 
                            - 6 * _shear_modulus_o * pow(_xi_0, 2) + 24 * _shear_modulus_o));
  Real denominator = 2 * (pow(_xi_0, 2) - 3);
  
  // Calculate gamma_r
  Real gamma_r = (term1 - term2_sqrt) / denominator;
  
  //save
  _gamma_damaged_r[_qp] = gamma_r;
}

void
DamageBreakageMaterial::computecoefficients()
{

  // Get xi_0 value
  const Real _xi_0 = _use_xi0_aux ? (*_xi0_aux)[_qp] : _xi0_value;
  // Get shear_modulus_o value
  const Real _shear_modulus_o = _use_shear_modulus_o_aux ? (*_shear_modulus_o_aux)[_qp] : _shear_modulus_o_value;

  //compute xi_1
  Real _xi_1 = _xi_0 + sqrt( pow(_xi_0 , 2) + 2 * _shear_modulus_o / _lambda_o );

  //compute alpha_cr | xi = 0
  Real alpha_cr_xi0 = alphacr_root1(0);

  //compute mu_cr
  Real mu_cr = _shear_modulus_o + alpha_cr_xi0 * _xi_0 * _gamma_damaged_r[_qp];

  //a0
  Real a0 = _chi * mu_cr;

  //a1
  Real numerator_a1 = -2 * _chi * mu_cr * pow(_xi_1, 3) + 6 * _chi * mu_cr * _xi_1 * pow(_xi_d, 2) - 4 * _chi * mu_cr * pow(_xi_d, 3)
                      - 2 * _gamma_damaged_r[_qp] * pow(_xi_1, 3) * _xi_d + 2 * _gamma_damaged_r[_qp] * pow(_xi_1, 3) * _xi_0
                      + _lambda_o * pow(_xi_1, 3) * pow(_xi_d, 2) + 2 * _shear_modulus_o * pow(_xi_1, 3);
  Real denominator_a1 = 2 * pow(_xi_1, 3) * _xi_d - 4 * pow(_xi_1, 2) * pow(_xi_d, 2) + 2 * _xi_1 * pow(_xi_d, 3);
  Real a1 = numerator_a1 / denominator_a1;

  //a2
  Real numerator_a2 = 2 * _chi * mu_cr * pow(_xi_1, 3) - 3 * _chi * mu_cr * pow(_xi_1, 2) * _xi_d + _chi * mu_cr * pow(_xi_d, 3)
                       + 2 * _gamma_damaged_r[_qp] * pow(_xi_1, 3) * _xi_d - 2 * _gamma_damaged_r[_qp] * pow(_xi_1, 3) * _xi_0
                       - _lambda_o * pow(_xi_1, 3) * pow(_xi_d, 2) - 2 * _shear_modulus_o * pow(_xi_1, 3);
  Real denominator_a2 = pow(_xi_1, 4) * _xi_d - 2 * pow(_xi_1, 3) * pow(_xi_d, 2) + pow(_xi_1, 2) * pow(_xi_d, 3); 
  Real a2 = numerator_a2 / denominator_a2; 

  //a3
  Real numerator_a3 = -2 * _chi * mu_cr * pow(_xi_1, 2) + 4 * _chi * mu_cr * _xi_1 * _xi_d - 2 * _chi * mu_cr * pow(_xi_d, 2)
                       - 2 * _gamma_damaged_r[_qp] * pow(_xi_1, 2) * _xi_d + 2 * _gamma_damaged_r[_qp] * pow(_xi_1, 2) * _xi_0
                       + _lambda_o * pow(_xi_1, 2) * pow(_xi_d, 2) + 2 * _shear_modulus_o * pow(_xi_1, 2);
  Real denominator_a3 = 2 * pow(_xi_1, 4) * _xi_d - 4 * pow(_xi_1, 3) * pow(_xi_d, 2) + 2 * pow(_xi_1, 2) * pow(_xi_d, 3);
  Real a3 = numerator_a3 / denominator_a3; 

  //save
  _a0[_qp] = a0;
  _a1[_qp] = a1;
  _a2[_qp] = a2;
  _a3[_qp] = a3;

}

// Function to compute initial damage with time dependent material
// This function is the same as "InitialDamageCycleSim2D"
void
DamageBreakageMaterial::computeinitialdamage2D() {

  /**
   * Input Parameters:
   * _build_param_peak_value: peak value of the initial damage
   * _build_param_sigma: standard deviation of the Gaussian distribution
   * _build_param_len_of_fault: length of the fault
   * _perturbation_build_param_use_damage_perturb: use damage perturbation flag
   */

  Real x_coord = _q_point[_qp](0); //along the strike direction
  Real y_coord = _q_point[_qp](1); //along the normal direction

  Real alpha_o = 0;

  Real r = 0.0;
  Real sigma = _build_param_sigma;
  if (x_coord > -0.5*_build_param_len_of_fault and x_coord < 0.5*_build_param_len_of_fault ){
    r = y_coord;
    alpha_o = std::max(_build_param_peak_value * std::exp(-1.0*(std::pow(r,2))/(sigma*sigma)),0.0);
  }
  else if (x_coord <= -0.5*_build_param_len_of_fault ){
    r = std::sqrt((y_coord - 0) * (y_coord - 0) + (x_coord - (-0.5*_build_param_len_of_fault )) * (x_coord - (-0.5*_build_param_len_of_fault )));
    alpha_o = std::max(_build_param_peak_value * std::exp(-1.0*(std::pow(r,2))/(sigma*sigma)),0.0);
  }
  else if (x_coord >= 0.5*_build_param_len_of_fault ){
    r = std::sqrt((y_coord - 0) * (y_coord - 0) + (x_coord - (0.5*_build_param_len_of_fault)) * (x_coord - (0.5*_build_param_len_of_fault )));
    alpha_o = std::max(_build_param_peak_value * std::exp(-1.0*(std::pow(r,2))/(sigma*sigma)),0.0);
  }

  if (_perturbation_build_param_use_damage_perturb)
  {
    alpha_o = alpha_o + _damage_perturbation[_qp];
  }

  _initial_damage_time_dependent_mat[_qp] = alpha_o; 

}

// Function for compute damage perturbation with time dependent material
// This function is the same as "DamagePerturbationSquare2D"
void
DamageBreakageMaterial::computedamageperturbation2D(){

  /**
   * Input Parameters:
   * _nucl_center: center of the nucleation zone
   * _length: length of the nucleation zone
   * _peak_damage: peak damage value
   * _sigma: standard deviation of the Gaussian distribution
   * _thickness: thickness of the nucleation zone
   * _duration: duration of the nucleation zone
   */

  //Get coordinates
  //here no rotation is applied yet
  Real xcoord = _q_point[_qp](0); //strike
  Real ycoord = _q_point[_qp](1); //normal

  //Get damage increments
  //Note: if the damage perturbation is used, the _dt within the nucleation phase (_t < _duration) should be used as constant
  //Otherwise there could be a problem (need to think carefully how to handle this better), a possible solution
  //is to find the function describing total damage perturbation: dd = k t, and if dd > 1, then dd = 1, and t > _duration, then dd = 0
  Real damage_inc = _perturbation_build_param_peak_value / (_perturbation_build_param_duration / _dt);

  //Assign initial damage perturbation
  Real dalpha = 0.0;
  if ( _t <= _perturbation_build_param_duration ){
    if ( (xcoord >= _perturbation_build_param_nucl_center[0] - _perturbation_build_param_length / 2.0) && (xcoord <= _perturbation_build_param_nucl_center[0] + _perturbation_build_param_length / 2.0) && (ycoord >= _perturbation_build_param_nucl_center[1] - _perturbation_build_param_thickness / 2.0) && (ycoord <= _perturbation_build_param_nucl_center[1] + _perturbation_build_param_thickness / 2.0) ){
      dalpha = _damage_perturbation_old[_qp] + damage_inc*std::exp(-(xcoord*xcoord)/(2.0*_perturbation_build_param_sigma*_perturbation_build_param_sigma));
    }
    else{
      dalpha = _damage_perturbation_old[_qp];
    }
  }
  else{
    dalpha = _damage_perturbation_old[_qp];
  }
  _damage_perturbation[_qp] = dalpha;  

}

// Function for alpha_func_root1
Real 
DamageBreakageMaterial::alphacr_root1(Real xi) {
    
    // Get xi_0 value
    const Real _xi_0 = _use_xi0_aux ? (*_xi0_aux)[_qp] : _xi0_value;
    // Get shear_modulus_o value
    const Real _shear_modulus_o = _use_shear_modulus_o_aux ? (*_shear_modulus_o_aux)[_qp] : _shear_modulus_o_value;

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
    Real denominator = 2 * _gamma_damaged_r[_qp] * (3 * pow(xi, 2) - 6 * xi * _xi_0 + 4 * _xi_0 * _xi_0 - 3);
    return (term1 - term2) / denominator;
}

// Function for alpha_func_root2
Real 
DamageBreakageMaterial::alphacr_root2(Real xi) {
    
    // Get xi_0 value
    const Real _xi_0 = _use_xi0_aux ? (*_xi0_aux)[_qp] : _xi0_value;
    // Get shear_modulus_o value
    const Real _shear_modulus_o = _use_shear_modulus_o_aux ? (*_shear_modulus_o_aux)[_qp] : _shear_modulus_o_value;

    return 2 * _shear_modulus_o / (_gamma_damaged_r[_qp] * (xi - 2 * _xi_0));
}

void
DamageBreakageMaterial::computePrincipalStrainAndOrientation(
    RealVectorValue & strain_in_crack_dir)
{
  // The rotation tensor is ordered such that directions for pre-existing cracks appear first
  // in the list of columns.  For example, if there is one existing crack, its direction is in the
  // first column in the rotation tensor.

  std::vector<Real> eigval(3, 0.0);
  RankTwoTensor eigvec;

  //the option of using either plastic or elastic strain rate
  //default is to use elastic strain rate
  if (_use_total_strain_rate)
    _total_lagrange_strain_old[_qp].symmetricEigenvaluesEigenvectors(eigval, eigvec);
  else
    _elastic_strain_old[_qp].symmetricEigenvaluesEigenvectors(eigval, eigvec);

  // If the elastic strain is beyond the cracking strain, save the eigen vectors as
  // the rotation tensor. Reverse their order so that the third principal strain
  // (most tensile) will correspond to the first crack.
  // _crack_rotation[_qp].fillColumn(0, eigvec.column(2));
  // _crack_rotation[_qp].fillColumn(1, eigvec.column(1));
  // _crack_rotation[_qp].fillColumn(2, eigvec.column(0));

  strain_in_crack_dir(0) = eigval[2];
  strain_in_crack_dir(1) = eigval[1];
  strain_in_crack_dir(2) = eigval[0];
}