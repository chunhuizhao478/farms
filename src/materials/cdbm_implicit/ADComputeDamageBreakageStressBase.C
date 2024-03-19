//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADComputeDamageBreakageStressBase.h"
#include "RankTwoTensor.h"
#include "SymmetricRankTwoTensor.h"

//template <typename RankTwoTensor>
InputParameters
ADComputeDamageBreakageStressBase::validParams()
{
  InputParameters params = Material::validParams();
  params.addParam<std::string>("base_name",
                               "Optional parameter that allows the user to define "
                               "multiple mechanics material systems on the same "
                               "block, i.e. for multiple phases");
  params.suppressParameter<bool>("use_displaced_mesh");
  params.addRequiredParam<Real>("lambda_o", "initial lambda value (first lame constant) [Pa]");
  params.addRequiredParam<Real>("shear_modulus_o", "initial shear modulus value (second lame constant) [Pa]");
  params.addParam<std::vector<MaterialPropertyName>>(
      "extra_stress_names",
      std::vector<MaterialPropertyName>(),
      "Material property names of rank two tensors to be added to the stress.");
  params.addRequiredParam<Real>("xi_0", "strain invariants ratio: onset of damage evolution");
  params.addRequiredParam<Real>( "gamma_damaged_r", "coefficient of damage solid modulus");
  params.addRequiredCoupledVar("initial_alpha", "initial distribution of alpha");
  return params;
}

//template <typename RankTwoTensor>
ADComputeDamageBreakageStressBase::ADComputeDamageBreakageStressBase(const InputParameters & parameters)
  : Material(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _mechanical_strain(getADMaterialProperty<RankTwoTensor>(_base_name + "mechanical_strain")),
    _stress(declareADProperty<RankTwoTensor>(_base_name + "stress")),
    _elastic_strain(declareADProperty<RankTwoTensor>(_base_name + "elastic_strain")),
    _extra_stresses(getParam<std::vector<MaterialPropertyName>>("extra_stress_names").size()),
    _alpha_damagedvar(declareADProperty<Real>("alpha_damagedvar")),
    _B(declareADProperty<Real>("B")),
    _xi(declareADProperty<Real>("xi")),
    _I1(declareADProperty<Real>("I1")),
    _I2(declareADProperty<Real>("I2")),
    _lambda(declareADProperty<Real>("lambda")),
    _shear_modulus(declareADProperty<Real>("shear_modulus")),
    _gamma_damaged(declareADProperty<Real>("gamma_damaged")),
    _eps_p(declareADProperty<RankTwoTensor>("eps_p")),
    _eps_e(declareADProperty<RankTwoTensor>("eps_e")),
    _eps_total(declareADProperty<RankTwoTensor>("eps_total")),
    _eps_total_init(declareADProperty<RankTwoTensor>("eps_total_init")),
    _sts_total(declareADProperty<RankTwoTensor>("sts_total")),
    _static_initial_stress_tensor(getADMaterialProperty<RankTwoTensor>("static_initial_stress_tensor")),
    _lambda_o(getParam<Real>("lambda_o")),
    _shear_modulus_o(getParam<Real>("shear_modulus_o"))
{
  if (getParam<bool>("use_displaced_mesh"))
    mooseError("The stress calculator needs to run on the undisplaced mesh.");

  const std::vector<MaterialPropertyName> extra_stress_names =
      getParam<std::vector<MaterialPropertyName>>("extra_stress_names");
  for (MooseIndex(_extra_stresses) i = 0; i < _extra_stresses.size(); ++i)
    _extra_stresses[i] = &getMaterialProperty<RankTwoTensor>(extra_stress_names[i]);
}

//template <typename RankTwoTensor>
void
ADComputeDamageBreakageStressBase::initQpStatefulProperties()
{
  _elastic_strain[_qp].zero();
  _stress[_qp].zero(); 

}

//template <typename RankTwoTensor>
void
ADComputeDamageBreakageStressBase::computeQpProperties()
{
  computeQpStress();

  // Add in extra stress
  for (MooseIndex(_extra_stresses) i = 0; i < _extra_stresses.size(); ++i)
    _stress[_qp] += (*_extra_stresses[i])[_qp];
}

//template class ADComputeDamageBreakageStressBaseTempl<RankTwoTensor>;
//template class ADComputeDamageBreakageStressBaseTempl<SymmetricRankTwoTensor>;