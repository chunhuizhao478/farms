//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ComputeDamageBreakageStress3DDynamicCDBM.h"
#include "NestedSolve.h"
#include "FEProblem.h"

registerMooseObject("farmsApp", ComputeDamageBreakageStress3DDynamicCDBM);

InputParameters
ComputeDamageBreakageStress3DDynamicCDBM::validParams()
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
  params.addRequiredParam<Real>(    "D", "D value");
  
  //variable parameters
  params.addRequiredCoupledVar("alpha_grad_x", "damage variable gradient component in x computed from subApp");
  params.addRequiredCoupledVar("alpha_grad_y", "damage variable gradient component in y computed from subApp");
  params.addRequiredCoupledVar("alpha_grad_z", "damage variable gradient component in z computed from subApp");

  //porous flow coupling
  params.addParam<bool>("porous_flow_coupling", false, "Flag to indicate if this material is coupled to a porous flow module");
  params.addParam<Real>("crack_surface_roughness_correction_factor", -1.0, "the factor to correct the crack opening based on crack surface roughness");
  params.addParam<Real>("length_scale", -1.0, "the length scale for nonlocal eqstrain");
  params.addParam<Real>("intrinsic_permeability", -1.0, "the intrisic permeability for the material, the intrisic permeability is used to take diagonal of the permeability tensor");
  return params;
}

ComputeDamageBreakageStress3DDynamicCDBM::ComputeDamageBreakageStress3DDynamicCDBM(const InputParameters & parameters)
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
    _crack_rotation_old(getMaterialPropertyOldByName<RankTwoTensor>("crack_rotation")),
    _alpha_grad_x(coupledValue("alpha_grad_x")),
    _alpha_grad_y(coupledValue("alpha_grad_y")),
    _alpha_grad_z(coupledValue("alpha_grad_z")),
    _D(getParam<Real>("D")),
    _initial_damage(getMaterialPropertyByName<Real>("initial_damage")),
    _initial_breakage(getMaterialPropertyByName<Real>("initial_breakage")),
    _shear_stress_perturbation(getMaterialPropertyOldByName<Real>("shear_stress_perturbation")),
    _Cd_constant(getParam<Real>("Cd_constant")),
    _C1(getParam<Real>("C_1")),
    _C2(getParam<Real>("C_2")),
    _beta_width(getParam<Real>("beta_width")),
    _CdCb_multiplier(getParam<Real>("CdCb_multiplier")),
    _CBH_constant(getParam<Real>("CBH_constant")),
    _dim(_mesh.dimension()),
    _porous_flow_coupling(getParam<bool>("porous_flow_coupling")),
    _crack_surface_roughness_correction_factor(getParam<Real>("crack_surface_roughness_correction_factor")),
    _length_scale(getParam<Real>("length_scale")),
    _intrinsic_permeability(getParam<Real>("intrinsic_permeability")),
    //## DEBUG ##
    _elasticity_tensor_name(_base_name + "elasticity_tensor"),
    _elasticity_tensor(getMaterialPropertyByName<RankFourTensor>(_elasticity_tensor_name))
{
  if (_porous_flow_coupling && ( _crack_surface_roughness_correction_factor < 0.0 || _length_scale < 0.0 || _intrinsic_permeability < 0.0) ){
    mooseError("The crack surface roughness correction factor, length scale and intrinsic permeability must be positive when using the porous flow coupling.");
  }
}

void
ComputeDamageBreakageStress3DDynamicCDBM::initialSetup()
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
ComputeDamageBreakageStress3DDynamicCDBM::initQpStatefulProperties()
{
  _elastic_strain[_qp].zero();
  _stress[_qp].zero();
  // _alpha_damagedvar[_qp] = _initial_damage[_qp];
  _crack_rotation[_qp] = RankTwoTensor::Identity(); 

}

void
ComputeDamageBreakageStress3DDynamicCDBM::computeQpStress()
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

  /* 
  compute alpha and B parameters
  */

  /* compute alpha */
  //compute forcing term
  Real alpha_forcingterm;
  if ( _xi_old[_qp] >= _xi_0 && _xi_old[_qp] <= _xi_max ){
    alpha_forcingterm = (1 - _B_old[_qp]) * ( _Cd_constant * _I2_old[_qp] * ( _xi_old[_qp] - _xi_0 ) );
  }
  else if ( _xi_old[_qp] < _xi_0 && _xi_old[_qp] >= _xi_min ){
    alpha_forcingterm = (1 - _B_old[_qp]) * ( _C1 * std::exp(_alpha_damagedvar_old[_qp]/_C2) * _I2_old[_qp] * ( _xi_old[_qp] - _xi_0 ) );
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
  if ( alpha_out < _initial_damage[_qp]){ alpha_out = _initial_damage[_qp]; }
  else{}

  _alpha_damagedvar[_qp] = alpha_out;

  /* compute B */
  Real C_B = _CdCb_multiplier * _Cd_constant;

  //compute xi_1
  Real _xi_1 = _xi_0 + sqrt( pow(_xi_0 , 2) + 2 * _shear_modulus_o / _lambda_o );

  //alphacr function
  Real alphacr;
  if ( _xi_old[_qp] < _xi_0 ){ alphacr = 1.0;} 
  else if ( _xi_old[_qp] > _xi_0 && _xi_old[_qp] <= _xi_1 ){ alphacr = alphacr_root1(_xi_old[_qp],gamma_damaged_r);}
  else if ( _xi_old[_qp] > _xi_1 && _xi_old[_qp] <= _xi_max ){ alphacr = alphacr_root2(_xi_old[_qp],gamma_damaged_r); }
  else{std::cout<<"xi: "<<_xi_old[_qp]<<std::endl;mooseError("xi exceeds the maximum allowable range!");}

  //compute forcing func
  Real Prob = 1.0 / ( std::exp( (alphacr - _alpha_damagedvar_old[_qp]) / _beta_width ) + 1.0 );
  Real B_forcingterm;
  if ( _xi_old[_qp] >= _xi_d && _xi_old[_qp] <= _xi_max ){
    B_forcingterm = 1.0 * C_B * Prob * (1-_B_old[_qp]) * _I2_old[_qp] * (_xi_old[_qp] - _xi_d); //could heal if xi < xi_0
  }
  else if ( _xi_old[_qp] < _xi_d && _xi_old[_qp] >= _xi_min ){
    B_forcingterm = 1.0 * _CBH_constant * _I2_old[_qp] * ( _xi_old[_qp] - _xi_d ); //close healing
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
  RankTwoTensor eps_e = _mechanical_strain[_qp] - eps_p;

  /* convert stress perturbation to strain perturbation */
  // Real shear_strain_perturbation = 0.0;
  // if (_shear_stress_perturbation[_qp] != 0){
  //   shear_strain_perturbation = _shear_stress_perturbation[_qp] / (2 * shear_modulus_out);
  //   eps_e(0,1) += shear_strain_perturbation;
  //   eps_e(1,0) += shear_strain_perturbation;
  // }  

  //const Real epsilon = 1e-16;
  Real I1 = eps_e(0,0) + eps_e(1,1) + eps_e(2,2);
  Real I2 = eps_e(0,0) * eps_e(0,0) + eps_e(1,1) * eps_e(1,1) + eps_e(2,2) * eps_e(2,2) + 2 * eps_e(0,1) * eps_e(0,1) + 2 * eps_e(0,2) * eps_e(0,2) + 2 * eps_e(1,2) * eps_e(1,2);
  Real xi = I1/std::sqrt(I2) + 1e-12; //catch the nan error in the initial solve

  //Represent sigma (solid(s) + granular(b))
  RankTwoTensor sigma_s;
  RankTwoTensor sigma_b;
  RankTwoTensor sigma_total;
  RankTwoTensor sigma_d;
  const auto I = RankTwoTensor::Identity();

  /* Compute effective elasticity tensor */
  // C_ijkl = lambda_eff * delta_ij * delta_kl + 2 * mu_eff * (delta_ik * delta_jl + delta_il * delta_jk);
  // sigma_ij = lambda_eff * delta_ij * I1 + 2 * mu_eff * eps_ij;
  // lambda_eff = (1 - B) * (lambda_out - gamma_damaged_out / xi) + B * (2 * a2 + a1 / xi + 3 * a3 * xi)
  // mu_eff = (1 - B) * 1/2 * (2 * shear_modulus_out - gamma_damaged_out * xi) + B * 1/2 * (2 * a0 + a1 * xi - a3 * std::pow(xi, 3))
  Real lambda_eff = (1 - B_out) * (lambda_out - gamma_damaged_out / xi) + B_out * (2 * a2 + a1 / xi + 3 * a3 * xi);
  Real mu_eff = (1 - B_out) * 0.5 * (2 * shear_modulus_out - gamma_damaged_out * xi) + B_out * 0.5 * (2 * a0 + a1 * xi - a3 * std::pow(xi, 3));
  _Eeff[_qp].fillFromInputVector({lambda_eff, mu_eff}, RankFourTensor::symmetric_isotropic);

  /* Compute stress */
  // sigma_s = (lambda_out - gamma_damaged_out / xi) * I1 * RankTwoTensor::Identity() + (2 * shear_modulus_out - gamma_damaged_out * xi) * eps_e;
  // sigma_b = (2 * a2 + a1 / xi + 3 * a3 * xi) * I1 * RankTwoTensor::Identity() + (2 * a0 + a1 * xi - a3 * std::pow(xi, 3)) * eps_e;
  sigma_total = (1 - B_out) * sigma_s + B_out * sigma_b;

  /* Compute Principal Strains and Rotation Matrix */
  RealVectorValue strain_in_crack_dir; //principal strains
  computeCrackStrainAndOrientation(strain_in_crack_dir);

  /* Update solid bulk compliance */
  updateSolidBulkCompliance();

  /* Update effective permeability */
  updatePermeabilityForCracking();  

  sigma_d = sigma_total - 0.3333 * (sigma_total(0,0) + sigma_total(1,1) + sigma_total(2,2)) * I;

  _eps_total[_qp] = eps_p + eps_e;
  _eps_p[_qp] = eps_p;
  _eps_e[_qp] = eps_e;
  _I1[_qp] = I1;
  _I2[_qp] = I2;
  _xi[_qp] = xi;
  _sigma_d[_qp] = sigma_d;

  // stress = C * e
  _stress[_qp] = _Eeff[_qp] * _mechanical_strain[_qp];

  // Assign value for elastic strain, which is equal to the mechanical strain
  _elastic_strain[_qp] = _mechanical_strain[_qp];

  // Compute dstress_dstrain
  _Jacobian_mult[_qp] = _Eeff[_qp];

  // // Rotate the stress state to the current configuration
  // _stress[_qp] = sigma_total;

  // // Assign value for elastic strain, which is equal to the mechanical strain
  // _elastic_strain[_qp] = eps_e;

  // // Compute tangent
  // RankFourTensor tangent;
  // computeQpTangentModulus(tangent,I1,I2,xi,eps_e);
  // _Jacobian_mult[_qp] = tangent;

  //Compute equivalent strain rate
  RankTwoTensor epsilon_rate = (eps_p - _eps_p_old[_qp])/_dt;
  Real epsilon_eq = sqrt(2/3*(epsilon_rate(0,0)*epsilon_rate(0,0)+epsilon_rate(1,1)*epsilon_rate(1,1)+epsilon_rate(2,2)*epsilon_rate(2,2)+2*epsilon_rate(0,1)*epsilon_rate(0,1)+2*epsilon_rate(0,2)*epsilon_rate(0,2)+2*epsilon_rate(1,2)*epsilon_rate(1,2)));
  _epsilon_eq[_qp] = epsilon_eq;

}

Real 
ComputeDamageBreakageStress3DDynamicCDBM::computegammar()
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
ComputeDamageBreakageStress3DDynamicCDBM::computecoefficients(Real gamma_damaged_r)
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
ComputeDamageBreakageStress3DDynamicCDBM::alphacr_root1(Real xi, Real gamma_damaged_r) {
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
ComputeDamageBreakageStress3DDynamicCDBM::alphacr_root2(Real xi, Real gamma_damaged_r) {
    return 2 * _shear_modulus_o / (gamma_damaged_r * (xi - 2 * _xi_0));
}

// void
// ComputeDamageBreakageStress3DDynamicCDBM::computeQpTangentModulus(RankFourTensor & tangent, 
//                                                       Real I1, 
//                                                       Real I2, 
//                                                       Real xi, 
//                                                       RankTwoTensor Ee)
// {

//   /*
//   compute gammar, breakage coefficients
//   */
//   Real gamma_damaged_r = computegammar();
//   std::vector<Real> avec = computecoefficients(gamma_damaged_r);
//   Real a0 = avec[0];
//   Real a1 = avec[1];
//   Real a2 = avec[2];
//   Real a3 = avec[3];

//   const Real adjusted_I2 = (I2 <= 1e-12) ? 1e-12 : I2;
//   const RankTwoTensor identity = RankTwoTensor::Identity();

//   // Precompute dxidE tensor
//   RankTwoTensor dxidE_tensor;
//   for (unsigned int k = 0; k < 3; ++k)
//     for (unsigned int l = 0; l < 3; ++l)
//       dxidE_tensor(k, l) = (identity(k, l) * adjusted_I2 - I1 * Ee(k, l)) / std::pow(adjusted_I2, 1.5);

//   const RankTwoTensor dxim1dE_tensor = dxidE_tensor * (-1.0 / (xi * xi));

//   // Compute terms for dSedE
//   const Real lambda_term = _lambda[_qp] - _gamma_damaged[_qp] / xi;
//   const Real shear_term = 2.0 * _shear_modulus[_qp] - _gamma_damaged[_qp] * xi;

//   RankFourTensor term_se1 = identity.outerProduct(-_gamma_damaged[_qp] * I1 * dxim1dE_tensor);
//   RankFourTensor term_se2 = identity.outerProduct(identity) * lambda_term;
//   RankFourTensor term_se3 = Ee.outerProduct(-_gamma_damaged[_qp] * dxidE_tensor);
//   RankFourTensor term_se4 = RankFourTensor(RankFourTensor::initIdentityFour) * shear_term;

//   RankFourTensor dSedE = term_se1 + term_se2 + term_se3 + term_se4;

//   // Compute terms for dSbdE
//   const Real coeff2_b = 2.0 * a2 + a1 / xi + 3.0 * a3 * xi;
//   const Real coeff4_b = 2.0 * a0 + a1 * xi - a3 * xi * xi * xi;

//   RankFourTensor term_b1 = identity.outerProduct((a1 * dxim1dE_tensor + 3 * a3 * dxidE_tensor) * I1);
//   RankFourTensor term_b2 = identity.outerProduct(identity) * coeff2_b;
//   RankFourTensor term_b3 = Ee.outerProduct(a1 * dxidE_tensor - a3 * 3 * xi * xi * dxidE_tensor);
//   RankFourTensor term_b4 = RankFourTensor(RankFourTensor::initIdentityFour) * coeff4_b;

//   RankFourTensor dSbdE = term_b1 + term_b2 + term_b3 + term_b4;

//   // Combine and assign tangent
//   tangent = (1.0 - _B[_qp]) * dSedE + _B[_qp] * dSbdE;  

// }

void
ComputeDamageBreakageStress3DDynamicCDBM::computeCrackStrainAndOrientation(
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
ComputeDamageBreakageStress3DDynamicCDBM::updateSolidBulkCompliance()
{

  // If porous flow coupling is not enabled, return
  if (!_porous_flow_coupling)
  return;

  // bulk modulus 
  // solid bulk moduls is a measure of a material's resistance to uniform compression
  // solid bulk modulus = (1/9) * sum_i(1,2,3) * sum_j(1,2,3) C_iijj
  Real solid_bulk_modulus = 0.0;
  RankFourTensor Eeff = _Eeff[_qp];
  for (int i = 0; i < 3; i++){
    for (int j = 0; j < 3; j++){
      solid_bulk_modulus += Eeff(i,i,j,j);
    }
  }

  solid_bulk_modulus = 1.0/9.0 * solid_bulk_modulus;

  //compute compliance
  _solid_bulk_compliance_damaged[_qp] = 1.0 / solid_bulk_modulus;
}

void
ComputeDamageBreakageStress3DDynamicCDBM::updatePermeabilityForCracking()
{

  // If porous flow coupling is not enabled, return
  if (!_porous_flow_coupling)
    return;

  // Get transformation matrix
  const RankTwoTensor & R = _crack_rotation[_qp];

  // The first column of the rotation tensor is the direction of the maximum principal strain, extract it as the normal to the crack plane
  RealVectorValue normal_to_crack_plane = R.column(0);

  // Normalize the normal to the crack plane
  RealVectorValue normalized_normal_to_crack_plane = normal_to_crack_plane / normal_to_crack_plane.norm();

  // Compute the crack opening
  Real xid = _alpha_damagedvar[_qp] -  0.5;
  HeavisideFunction(xid);

  // Compute open crack width w_c
  Real w_c = std::norm(_length_scale * ( 1 + normalized_normal_to_crack_plane * _elastic_strain[_qp] * normalized_normal_to_crack_plane ) );

  // Compute closed crack width w_r
  Real w_r = std::sqrt(12 * ( 10 * _intrinsic_permeability) );

  // Compute the aperture w_h, take the maximum between w_c and w_r
  Real w_h = std::max(_crack_surface_roughness_correction_factor * w_c * xid, _crack_surface_roughness_correction_factor * w_r * xid);

  //Compute the effective permeability
  RankTwoTensor perm_frac = (w_h * w_h / 12.0) * (RankTwoTensor::Identity() - RankTwoTensor::outerProduct(normalized_normal_to_crack_plane, normalized_normal_to_crack_plane));

  //Compute the intrinsic permeability
  RankTwoTensor perm_intrinsic = _intrinsic_permeability * RankTwoTensor::Identity();

  //Compute the effective permeability
  _effective_perm[_qp] = perm_intrinsic + _alpha_damagedvar[_qp] * perm_frac;

  // Ensure effective permeability is not less than intrinsic permeability
  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
  {
      for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
      {
          if (_effective_perm[_qp](i,j) < perm_intrinsic(i,j))
          {
            _effective_perm[_qp](i,j) = perm_intrinsic(i,j);
          }
      }
  }

}

void
ComputeDamageBreakageStress3DDynamicCDBM::HeavisideFunction(Real & x)
{
  if (x > 0.0)
    x = 1.0;
  else
    x = 0.0;
}