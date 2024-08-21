//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ComputeDamageBreakageStress3D.h"
#include "NestedSolve.h"
#include "FEProblem.h"

registerMooseObject("farmsApp", ComputeDamageBreakageStress3D);

InputParameters
ComputeDamageBreakageStress3D::validParams()
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
  params.addRequiredParam<Real>(            "xi_1", "critical point of three phases");
  params.addRequiredParam<Real>(          "xi_min", "strain invariants ratio: minimum allowable value");
  params.addRequiredParam<Real>(          "xi_max", "strain invariants ratio: maximum allowable value");
  params.addRequiredParam<Real>(              "a0", "parameters in granular states");
  params.addRequiredParam<Real>(              "a1", "parameters in granular states");
  params.addRequiredParam<Real>(              "a2", "parameters in granular states");
  params.addRequiredParam<Real>(              "a3", "parameters in granular states");
  params.addRequiredParam<Real>( "gamma_damaged_r", "coefficient of damage solid modulus");
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

  return params;
}

ComputeDamageBreakageStress3D::ComputeDamageBreakageStress3D(const InputParameters & parameters)
  : ComputeDamageBreakageStressBase3D(parameters),
    _xi_0(getParam<Real>("xi_0")),
    _xi_d(getParam<Real>("xi_d")),
    _xi_1(getParam<Real>("xi_1")),
    _xi_min(getParam<Real>("xi_min")),
    _xi_max(getParam<Real>("xi_max")),
    _a0(getParam<Real>("a0")),
    _a1(getParam<Real>("a1")),
    _a2(getParam<Real>("a2")),
    _a3(getParam<Real>("a3")),
    _gamma_damaged_r(getParam<Real>("gamma_damaged_r")),
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
    _alpha_grad_x(coupledValue("alpha_grad_x")),
    _alpha_grad_y(coupledValue("alpha_grad_y")),
    _alpha_grad_z(coupledValue("alpha_grad_z")),
    _density_old(getMaterialPropertyOldByName<Real>("density")),
    _D(getParam<Real>("D")),
    _static_initial_stress_tensor(getMaterialPropertyByName<RankTwoTensor>("static_initial_stress_tensor")),
    _static_initial_strain_tensor(getMaterialPropertyByName<RankTwoTensor>("static_initial_strain_tensor")),
    _I1_initial(getMaterialPropertyByName<Real>("I1_initial")),
    _I2_initial(getMaterialPropertyByName<Real>("I2_initial")),
    _xi_initial(getMaterialPropertyByName<Real>("xi_initial")),
    _initial_damage(getMaterialPropertyByName<Real>("initial_damage")),
    _Cd_constant(getParam<Real>("Cd_constant")),
    _C1(getParam<Real>("C_1")),
    _C2(getParam<Real>("C_2")),
    _beta_width(getParam<Real>("beta_width")),
    _CdCb_multiplier(getParam<Real>("CdCb_multiplier")),
    _CBH_constant(getParam<Real>("CBH_constant"))
{
}

void
ComputeDamageBreakageStress3D::initialSetup()
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
ComputeDamageBreakageStress3D::initQpStatefulProperties()
{
  _elastic_strain[_qp].zero();
  _stress[_qp].zero();

}

void
ComputeDamageBreakageStress3D::computeQpStress()
{ 
  
  if ( _t == 0 ){

    _alpha_damagedvar[_qp] = _initial_damage[_qp];
    _B[_qp] = 0.0;

    _xi[_qp] = _xi_initial[_qp];
    _I1[_qp] = _I1_initial[_qp];
    _I2[_qp] = _I2_initial[_qp];

    /// lambda (first lame const)
    _lambda[_qp] = _lambda_o;
    /// mu (shear modulus)
    _shear_modulus[_qp] = _shear_modulus_o + _initial_damage[_qp] * _xi_0 * _gamma_damaged_r;
    /// gamma_damaged (damage modulus)
    _gamma_damaged[_qp] = _initial_damage[_qp] * _gamma_damaged_r;
    /// initialize strain tensor eps_p eps_e eps_total xi I1 I2
    _eps_p[_qp].zero();
    _eps_e[_qp] = _static_initial_strain_tensor[_qp];
    _eps_total[_qp] = _static_initial_strain_tensor[_qp];
    _sts_total[_qp] = _static_initial_stress_tensor[_qp];
  
  }
  else{
  
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
  if ( alpha_out < _initial_damage[_qp] ){ alpha_out = _initial_damage[_qp]; }
  else{}

  _alpha_damagedvar[_qp] = alpha_out;

  //grad_alpha
  // Real alpha_grad_x = _alpha_grad_x[_qp];
  // Real alpha_grad_y = _alpha_grad_y[_qp];
  // Real alpha_grad_z = _alpha_grad_z[_qp];
  // Real D = _D;

  /* compute B */
  Real C_B = _CdCb_multiplier * _Cd_constant;

  //alphacr function
  Real alphacr;
  if ( _xi_old[_qp] < _xi_0 ){ alphacr = 1.0;} 
  else if ( _xi_old[_qp] > _xi_0 && _xi_old[_qp] <= _xi_1 ){ alphacr = ((_xi_old[_qp]*2.76e5-7.100521107637101e2*_xi_old[_qp]*7.5e2-7.100521107637101e2*1.4e3-7.100521107637101e2*std::pow(_xi_old[_qp],3)*1.25e2+std::pow(_xi_old[_qp],3)*4.6e4+std::sqrt((7.100521107637101e2*3.68e2-3.19799e5)*(_xi_old[_qp]*-1.44e3-std::pow(_xi_old[_qp],2)*2.1e3+std::pow(_xi_old[_qp],3)*5.6e2+std::pow(_xi_old[_qp],4)*3.0e2+std::pow(_xi_old[_qp],6)*2.5e1+3.576e3)*(-3.590922148807814e-1))*5.9e1+5.152e5)*(-5.9e1/4.0))/(_xi_old[_qp]*3.837588e7-7.100521107637101e2*_xi_old[_qp]*4.416e4+7.100521107637101e2*4.048e3-7.100521107637101e2*std::pow(_xi_old[_qp],2)*2.76e4+std::pow(_xi_old[_qp],2)*2.3984925e7-3.517789e6);}
  else if ( _xi_old[_qp] > _xi_1 && _xi_old[_qp] <= _xi_max ){ alphacr = 6.0e10/(7.100521107637101e2*1.627118644067797e8+_xi_old[_qp]*(7.100521107637101e2*1.016949152542373e8-3.742372881355932e10)-5.987796610169492e10);}
  else{std::cout<<"xi: "<<_xi_old[_qp]<<std::endl;mooseError("xi exceeds the maximum allowable range!");}

  //compute forcing func
  Real Prob = 1.0 / ( std::exp( (alphacr - _alpha_damagedvar_old[_qp]) / _beta_width ) + 1.0 );
  Real B_forcingterm;
  // if ( _xi_old[_qp] >= _xi_d && _xi_old[_qp] <= _xi_max ){
  //   B_forcingterm = 1.0 * C_B * Prob * (1-_B_old[_qp]) * _I2_old[_qp] * (_xi_old[_qp] - _xi_d); //could heal if xi < xi_0
  // }
  // else if ( _xi_old[_qp] < _xi_d && _xi_old[_qp] >= _xi_min ){
  //   B_forcingterm = 1.0 * _CBH_constant * _I2_old[_qp] * ( _xi_old[_qp] - _xi_d ); //close healing
  // }
  // else{
  //   mooseError("xi_old is OUT-OF-RANGE!.");
  // }

  //ggw183
  if ( _xi_old[_qp] >= _xi_d && _xi_old[_qp] <= _xi_max ){
    B_forcingterm = 1.0 * C_B * Prob * (1-_B_old[_qp]) * _I2_old[_qp] * ((_shear_modulus_old[_qp]-_a0)-(_a1+_gamma_damaged_old[_qp])*_xi_old[_qp]+(0.5*_lambda_o-_a2)*_xi_old[_qp]*_xi_old[_qp]-(_a3)*_xi_old[_qp]*_xi_old[_qp]*_xi_old[_qp])/(1000); //could heal if xi < xi_0
  }
  else if ( _xi_old[_qp] < _xi_d && _xi_old[_qp] >= _xi_min ){
    B_forcingterm = 1.0 * _CBH_constant * _I2_old[_qp] * ((_shear_modulus_old[_qp]-_a0)-(_a1+_gamma_damaged_old[_qp])*_xi_old[_qp]+(0.5*_lambda_o-_a2)*_xi_old[_qp]*_xi_old[_qp]-(_a3)*_xi_old[_qp]*_xi_old[_qp]*_xi_old[_qp])/(1000);
  }
  else{
    mooseError("xi_old is OUT-OF-RANGE!.");
  }

  Real B_out = _B_old[_qp] + _dt * B_forcingterm;

  //check breakage within range
  if ( B_out < 0 ){ B_out = 0.0; }
  else if ( B_out > 1 ){ B_out = 1.0; }
  else{}   

  //save alpha and B
  _alpha_damagedvar[_qp] = alpha_out;
  _B[_qp] = B_out;

  /*
    update modulus
  */

  //lambda, shear_modulus, gamma_damaged are updated
  Real lambda_out = _lambda_o;
  Real shear_modulus_out = _shear_modulus_o + alpha_out * _xi_0 * _gamma_damaged_r;
  Real gamma_damaged_out = alpha_out * _gamma_damaged_r;

  //save
  _lambda[_qp] = lambda_out;
  _shear_modulus[_qp] = shear_modulus_out;
  _gamma_damaged[_qp] = gamma_damaged_out;

  /*
    compute strain
  */

  /* compute strain */
  RankTwoTensor eps_p = _eps_p_old[_qp] + _dt * _C_g * std::pow(_B_old[_qp],_m1) * _sigma_d_old[_qp];
  RankTwoTensor eps_e = _mechanical_strain[_qp] - eps_p;
  Real I1 = eps_e(0,0) + eps_e(1,1) + eps_e(2,2);
  Real I2 = eps_e(0,0) * eps_e(0,0) + eps_e(1,1) * eps_e(1,1) + eps_e(2,2) * eps_e(2,2) + 2 * eps_e(0,1) * eps_e(0,1) + 2 * eps_e(0,2) * eps_e(0,2) + 2 * eps_e(1,2) * eps_e(1,2);
  Real xi = I1/std::sqrt(I2);

  //Represent sigma (solid(s) + granular(b))
  RankTwoTensor sigma_s;
  RankTwoTensor sigma_b;
  RankTwoTensor sigma_total;
  RankTwoTensor sigma_d;
  const auto I = RankTwoTensor::Identity();

  sigma_s(0,0) = ( _lambda_o - gamma_damaged_out / xi ) * I1 + ( 2 * shear_modulus_out - gamma_damaged_out * xi ) * eps_e(0,0);
  sigma_b(0,0) = ( 2 * _a2 + _a1 / xi + 3 * _a3 * xi ) * I1 + ( 2 * _a0 + _a1 * xi - _a3 * std::pow(xi,3) ) * eps_e(0,0);
  sigma_s(1,1) = ( _lambda_o - gamma_damaged_out / xi ) * I1 + ( 2 * shear_modulus_out - gamma_damaged_out * xi ) * eps_e(1,1);
  sigma_b(1,1) = ( 2 * _a2 + _a1 / xi + 3 * _a3 * xi ) * I1 + ( 2 * _a0 + _a1 * xi - _a3 * std::pow(xi,3) ) * eps_e(1,1);
  sigma_s(2,2) = ( _lambda_o - gamma_damaged_out / xi ) * I1 + ( 2 * shear_modulus_out - gamma_damaged_out * xi ) * eps_e(2,2);
  sigma_b(2,2) = ( 2 * _a2 + _a1 / xi + 3 * _a3 * xi ) * I1 + ( 2 * _a0 + _a1 * xi - _a3 * std::pow(xi,3) ) * eps_e(2,2);
  sigma_s(0,1)  = ( 2 * shear_modulus_out - gamma_damaged_out * xi    ) * eps_e(0,1); sigma_s(1,0)  = ( 2 * shear_modulus_out - gamma_damaged_out * xi    ) * eps_e(1,0);
  sigma_b(0,1)  = ( 2 * _a0 + _a1 * xi - _a3 * std::pow(xi,3) ) * eps_e(0,1); sigma_b(1,0)  = ( 2 * _a0 + _a1 * xi - _a3 * std::pow(xi,3) ) * eps_e(1,0);
  sigma_s(0,2)  = ( 2 * shear_modulus_out - gamma_damaged_out * xi    ) * eps_e(0,2); sigma_s(2,0)  = ( 2 * _a0 + _a1 * xi - _a3 * std::pow(xi,3) ) * eps_e(2,0);
  sigma_b(0,2)  = ( 2 * _a0 + _a1 * xi - _a3 * std::pow(xi,3) ) * eps_e(0,2); sigma_b(2,0)  = ( 2 * _a0 + _a1 * xi - _a3 * std::pow(xi,3) ) * eps_e(2,0);
  sigma_s(1,2)  = ( 2 * shear_modulus_out - gamma_damaged_out * xi    ) * eps_e(1,2); sigma_s(2,1)  = ( 2 * _a0 + _a1 * xi - _a3 * std::pow(xi,3) ) * eps_e(2,1);
  sigma_b(1,2)  = ( 2 * _a0 + _a1 * xi - _a3 * std::pow(xi,3) ) * eps_e(1,2); sigma_b(2,1)  = ( 2 * _a0 + _a1 * xi - _a3 * std::pow(xi,3) ) * eps_e(2,1);
  
  sigma_total = (1 - _B_old[_qp]) * sigma_s + _B_old[_qp] * sigma_b;
  sigma_d = sigma_total - 1/3 * (sigma_total(0,0) + sigma_total(1,1) + sigma_total(2,2)) * I;

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
  _elastic_strain[_qp] = _mechanical_strain[_qp];
  }
}

///Function: deltaij
Real 
ComputeDamageBreakageStress3D::deltaij(int i, int j)
{
  Real deltaij_out = 0.0;
  if ( i == j )
  {
    deltaij_out = 1.0;
  }
  else
  {
    deltaij_out = 0.0;
  }
  return deltaij_out;
}

/// Function: epsilonij - take component of elastic strain
Real 
ComputeDamageBreakageStress3D::epsilonij(int i, 
                                       int j,
                                       Real eps11e_in,
                                       Real eps22e_in,
                                       Real eps12e_in,
                                       Real eps33e_in,
                                       Real eps13e_in,
                                       Real eps23e_in)
{ 
  //take elastic strain
  Real eps_e_out = 0.0;
  if ( i == 1 && j == 1 )
  {
    eps_e_out = eps11e_in;
  }
  else if ( i == 2 && j == 2 )
  {
    eps_e_out = eps22e_in;
  }
  else if ( i == 3 && j == 3 )
  {
    eps_e_out = eps33e_in;
  }
  else if ( i == 1 && j == 2 )
  {
    eps_e_out = eps12e_in;
  }
  else if ( i == 1 && j == 3 )
  {
    eps_e_out = eps13e_in;
  }
  else if ( i == 2 && j == 3 )
  {
    eps_e_out = eps23e_in;
  }
  return eps_e_out;
}

/// Function: grad_alpha
//Note: no grad in 33 direction
Real 
ComputeDamageBreakageStress3D::grad_alpha(int i, 
                                          Real alpha_grad_x,
                                          Real alpha_grad_y,
                                          Real alpha_grad_z)
{ 
  //take alpha grad
  Real alpha_grad_out = 0.0;
  if ( i == 1 )
  {
    alpha_grad_out = alpha_grad_x;
  }
  else if ( i == 2 )
  {
    alpha_grad_out = alpha_grad_y;
  }
  else if ( i == 3 )
  {
    alpha_grad_out = alpha_grad_z;
  }
  return alpha_grad_out;
}

/// Function: compute stress components
Real 
ComputeDamageBreakageStress3D::computeStressComps(int i, 
                                                int j,
                                                Real xi_in,
                                                Real I1_in,
                                                Real B_in,
                                                Real lambda_in,
                                                Real gamma_damaged_in,
                                                Real shear_modulus_in,
                                                Real eps11e_in,
                                                Real eps22e_in,
                                                Real eps12e_in,
                                                Real eps33e_in,
                                                Real eps13e_in,
                                                Real eps23e_in,
                                                Real alpha_grad_x,
                                                Real alpha_grad_y,
                                                Real alpha_grad_z,
                                                Real D)
{
  //Retrieve parameters
  Real xi = xi_in;
  Real I1 = I1_in;
  Real B = B_in;
  Real lambda = lambda_in;
  Real gamma_damaged = gamma_damaged_in;
  Real shear_modulus = shear_modulus_in;

  //stress comps
  Real stresscomp_s_out;
  Real stresscomp_b_out;
  Real stresscomp;

  stresscomp_s_out = (lambda - gamma_damaged / xi) * I1 * deltaij(i,j) + ( 2 * shear_modulus - gamma_damaged * xi ) * epsilonij(i,j,eps11e_in,eps22e_in,eps12e_in,eps33e_in,eps13e_in,eps23e_in) - D * grad_alpha(i,alpha_grad_x,alpha_grad_y,alpha_grad_z) * grad_alpha(j,alpha_grad_x,alpha_grad_y,alpha_grad_z);
  stresscomp_b_out = (2 * _a2 + _a1 / xi + 3 * _a3 * xi) * I1 * deltaij(i,j) + ( 2 * _a0 + _a1 * xi - _a3 * pow(xi,3) ) * epsilonij(i,j,eps11e_in,eps22e_in,eps12e_in,eps33e_in,eps13e_in,eps23e_in);

  //
  stresscomp = ( 1 - B ) * stresscomp_s_out + B * stresscomp_b_out;

  return stresscomp;

}

///Function to compute initial strain based on initial stress 
void
ComputeDamageBreakageStress3D::setupInitial()
{
  // ///For isotropic material with all components of stress subject to small strain we consider 
  // ///tensile/compressive stress leads to only tensile/compressive strain, shear stress produce
  // ///shear strain: 
  // ///eps_ii = 1/E * ( sigma_ii - nu * ( sigma_jj + sigma_kk) )
  // ///eps_ij = 1/G * sigma_ij

  // /// lambda (first lame const)
  // _lambda[_qp] = _lambda_o;
  // /// mu (shear modulus)
  // _shear_modulus[_qp] = _shear_modulus_o;
  // /// gamma_damaged (damage modulus)
  // _gamma_damaged[_qp] = 0.0;

  // //allpha, B
  // _alpha_damagedvar[_qp] = 0.0;
  // _B[_qp] = 0.0;

  // //Convert (lambda_o,shear_modulus_o) to (youngs_modulus_o,poisson_ratio_o)
  // Real youngs_modulus_o = _shear_modulus_o * ( 3 * _lambda_o + 2 * _shear_modulus_o ) / ( _lambda_o + _shear_modulus_o );
  // Real poisson_ratio_o = _lambda_o / ( 2 * ( _lambda_o + _shear_modulus_o ));

  // //Convert (lambda_o,shear_modulus_o) to (shear_wave_speed_o,pressure_wave_speed_o)
  // // Real density_o = _density_old[_qp];
  // // Real shear_wave_speed_o = sqrt( ( _shear_modulus_o ) / ( density_o ) );
  // // Real pressure_wave_speed_o = sqrt( ( _lambda_o + 2 * _shear_modulus_o ) / ( density_o ) );

  // //save
  // // _shear_wave_speed[_qp] = shear_wave_speed_o;
  // // _pressure_wave_speed[_qp] = pressure_wave_speed_o;

  // //Get stress components
  // RankTwoTensor stress_initial = _static_initial_stress_tensor[_qp];

  // //Get strain components
  // RankTwoTensor strain_initial = _static_initial_strain_tensor[_qp];

  // //Compute strain components
  // Real eps11_init = strain_initial(0,0);
  // Real eps22_init = strain_initial(0,1);
  // Real eps12_init = strain_initial(0,2);
  // Real eps13_init = strain_initial(1,2);
  // Real eps23_init = strain_initial(1,1);
  // Real eps33_init = strain_initial(2,2); 

  // //Compute xi, I1, I2
  // Real I1_init = eps11_init + eps22_init + eps33_init;
  // Real I2_init = eps11_init * eps11_init + eps22_init * eps22_init + eps33_init * eps33_init + 2 * eps12_init * eps12_init + 2 * eps13_init * eps13_init + 2 * eps23_init * eps23_init;
  // Real xi_init = I1_init / sqrt( I2_init );

  // //Compute eps
  // //eps_p
  // _eps_p[_qp](0,0) = 0.0; _eps_p[_qp](0,1) = 0.0; _eps_p[_qp](0,2) = 0.0;
  // _eps_p[_qp](1,0) = 0.0; _eps_p[_qp](1,1) = 0.0; _eps_p[_qp](1,2) = 0.0;
  // _eps_p[_qp](2,0) = 0.0; _eps_p[_qp](2,1) = 0.0; _eps_p[_qp](2,2) = 0.0;
  // //eps_e
  // _eps_e[_qp](0,0) = eps11_init; _eps_e[_qp](0,1) = eps12_init; _eps_e[_qp](0,2) = eps13_init;
  // _eps_e[_qp](1,0) = eps12_init; _eps_e[_qp](1,1) = eps22_init; _eps_e[_qp](1,2) = eps23_init;
  // _eps_e[_qp](2,0) = eps13_init; _eps_e[_qp](2,1) = eps23_init; _eps_e[_qp](2,2) = eps33_init;
  // //eps_total
  // _eps_total[_qp](0,0) = eps11_init; _eps_total[_qp](0,1) = eps12_init; _eps_total[_qp](0,2) = eps13_init;
  // _eps_total[_qp](1,0) = eps12_init; _eps_total[_qp](1,1) = eps22_init; _eps_total[_qp](1,2) = eps23_init;
  // _eps_total[_qp](2,0) = eps13_init; _eps_total[_qp](2,1) = eps23_init; _eps_total[_qp](2,2) = eps33_init;
  // //sts_total
  // _sts_total[_qp] = stress_initial;

  // //I1
  // _I1[_qp] = I1_init;
  // //I2
  // _I2[_qp] = I2_init;
  // //xi
  // _xi[_qp] = xi_init;
}