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
    _sts_total(declareProperty<RankTwoTensor>("sts_total")),
    _sigma_d(declareProperty<RankTwoTensor>("sigma_d")),
    _epsilon_eq(declareProperty<Real>("epsilon_eq")),
    _lambda_o(getParam<Real>("lambda_o")),
    _shear_modulus_o(getParam<Real>("shear_modulus_o"))
{
}

void
ComputeGeneralDamageBreakageStressBase3D::computeQpProperties()
{
  computeQpStress();
}
