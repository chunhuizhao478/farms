//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ComputeGeneralDamageBreakageStressBase3D.h"
#include "ComputeElasticityTensorBase.h"
#include "Function.h"

InputParameters
ComputeGeneralDamageBreakageStressBase3D::validParams()
{
  InputParameters params = Material::validParams();
  params.addRequiredParam<Real>("lambda_o", "initial lambda value (first lame constant) [Pa]");
  params.addRequiredParam<Real>("shear_modulus_o", "initial shear modulus value (second lame constant) [Pa]");
  params.addParam<std::string>("base_name",
                               "Optional parameter that allows the user to define "
                               "multiple mechanics material systems on the same "
                               "block, i.e. for multiple phases");
  return params;
}

ComputeGeneralDamageBreakageStressBase3D::ComputeGeneralDamageBreakageStressBase3D(const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _mechanical_strain(getMaterialPropertyByName<RankTwoTensor>(_base_name + "mechanical_strain")),
    _stress(declareProperty<RankTwoTensor>(_base_name + "stress")),
    _elastic_strain(declareProperty<RankTwoTensor>(_base_name + "elastic_strain")),
    // _extra_stress(getDefaultMaterialProperty<RankTwoTensor>(_base_name + "extra_stress")),
    _Jacobian_mult(declareProperty<RankFourTensor>(_base_name + "Jacobian_mult")),
    _alpha_damagedvar(declareProperty<Real>("alpha_damagedvar")),
    _B(declareProperty<Real>("B")),
    _xi(declareProperty<Real>("xi")),
    _I1(declareProperty<Real>("I1")),
    _I2(declareProperty<Real>("I2")),
    _lambda(declareProperty<Real>("lambda")),
    _shear_modulus(declareProperty<Real>("shear_modulus")),
    _gamma_damaged(declareProperty<Real>("gamma_damaged")),
    _eps_p(declareProperty<RankTwoTensor>("eps_p")),
    _eps_e(declareProperty<RankTwoTensor>("eps_e")),
    _eps_total(declareProperty<RankTwoTensor>("eps_total")),
    // _eps_total_init(declareProperty<RankTwoTensor>("eps_total_init")),
    _sts_total(declareProperty<RankTwoTensor>("sts_total")),
    // _shear_wave_speed(declareProperty<Real>("shear_wave_speed")),
    // _pressure_wave_speed(declareProperty<Real>("pressure_wave_speed")),
    _static_initial_stress_tensor(getMaterialPropertyByName<RankTwoTensor>("static_initial_stress_tensor")),
    // _principal_strain(declareProperty<Real>("principal_strain")),
    // _eqv_plastic_strain(declareProperty<Real>("eqv_plastic_strain")),
    // _plastic_work(declareProperty<Real>("plastic_work")),
    // _plastic_work_rate(declareProperty<Real>("plastic_work_rate")),
    // _resolved_sigma_11(declareProperty<Real>("resolved_sigma_11")),
    // _resolved_sigma_12(declareProperty<Real>("resolved_sigma_12")),
    // _resolved_sigma_22(declareProperty<Real>("resolved_sigma_22")),
    // _resolved_sigma_33(declareProperty<Real>("resolved_sigma_33")),
    // _resolved_epsp_rate_11(declareProperty<Real>("resolved_epsp_rate_11")),
    // _resolved_epsp_rate_12(declareProperty<Real>("resolved_epsp_rate_12")),
    // _resolved_epsp_rate_22(declareProperty<Real>("resolved_epsp_rate_22")),
    // _resolved_epsp_rate_33(declareProperty<Real>("resolved_epsp_rate_33")),
    _lambda_o(getParam<Real>("lambda_o")),
    _shear_modulus_o(getParam<Real>("shear_modulus_o"))
{
}

void
ComputeGeneralDamageBreakageStressBase3D::initQpStatefulProperties()
{
  _elastic_strain[_qp].zero();
  _stress[_qp].zero();

  // /// lambda (first lame const)
  // _lambda[_qp] = _lambda_o;
  // /// mu (shear modulus)
  // _shear_modulus[_qp] = _shear_modulus_o;
  // /// gamma_damaged (damage modulus)
  // _gamma_damaged[_qp] = 0.0;

  /// initialize strain tensor eps_p eps_e eps_total xi I1 I2
  ///computeInitialStrain();

}

void
ComputeGeneralDamageBreakageStressBase3D::computeQpProperties()
{
  computeQpStress();

  // Add in extra stress
  // _stress[_qp] += _extra_stress[_qp];
}

// ///Function to compute initial strain based on initial stress 
// void
// ComputeGeneralDamageBreakageStressBase3D::computeInitialStrain()
// {
//   ///For isotropic material with all components of stress subject to small strain we consider 
//   ///tensile/compressive stress leads to only tensile/compressive strain, shear stress produce
//   ///shear strain: 
//   ///eps_ii = 1/E * ( sigma_ii - nu * ( sigma_jj + sigma_kk) )
//   ///eps_ij = 1/G * sigma_ij

//   //Convert (lambda_o,shear_modulus_o) to (youngs_modulus_o,poisson_ratio_o)
//   Real youngs_modulus_o = _shear_modulus_o * ( 3 * _lambda_o + 2 * _shear_modulus_o ) / ( _lambda_o + _shear_modulus_o );
//   Real poisson_ratio_o = _lambda_o / ( 2 * ( _lambda_o + _shear_modulus_o ));

//   //Get stress components
//   RankTwoTensor stress_initial = _static_initial_stress_tensor[_qp];
//   Real sts11_init = stress_initial(0,0);
//   Real sts12_init = stress_initial(0,1);
//   Real sts22_init = stress_initial(1,1);

//   //Compute strain components
//   Real eps11_init = 1.0 / youngs_modulus_o * ( sts11_init - poisson_ratio_o * ( sts22_init ) );
//   Real eps22_init = 1.0 / youngs_modulus_o * ( sts22_init - poisson_ratio_o * ( sts11_init ) );
//   Real eps12_init = 1.0 / _shear_modulus_o * sts12_init;

//   //Compute xi, I1, I2
//   Real I1_init = eps11_init + eps22_init;
//   Real I2_init = eps11_init * eps11_init + eps22_init * eps22_init + 2 * eps12_init * eps12_init;
//   Real xi_init = I1_init / sqrt( I2_init );

//   //Compute eps
//   //eps_p
//   _eps_p[_qp](0,0) = 0.0; _eps_p[_qp](0,1) = 0.0; 
//   _eps_p[_qp](1,0) = 0.0; _eps_p[_qp](1,1) = 0.0;
//   //eps_e
//   _eps_e[_qp](0,0) = eps11_init; _eps_e[_qp](0,1) = eps12_init; 
//   _eps_e[_qp](1,0) = eps12_init; _eps_e[_qp](1,1) = eps22_init;
//   //eps_total
//   _eps_total[_qp](0,0) = eps11_init; _eps_total[_qp](0,1) = eps12_init;
//   _eps_total[_qp](1,0) = eps12_init; _eps_total[_qp](1,1) = eps22_init;
//   //eps_total_init
//   _eps_total_init[_qp](0,0) = eps11_init; _eps_total_init[_qp](0,1) = eps12_init;
//   _eps_total_init[_qp](1,0) = eps12_init; _eps_total_init[_qp](1,1) = eps22_init;
//   //I1
//   _I1[_qp] = I1_init;
//   //I2
//   _I2[_qp] = I2_init;
//   //xi
//   _xi[_qp] = xi_init;
// }