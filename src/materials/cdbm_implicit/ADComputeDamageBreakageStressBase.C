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
    _B_breakagevar(declareADProperty<Real>("B_breakagevar")),
    _xi(declareADProperty<Real>("xi")),
    _I1(declareADProperty<Real>("I1")),
    _I2(declareADProperty<Real>("I2")),
    _shear_modulus(declareADProperty<Real>("shear_modulus")),
    _gamma_damaged(declareADProperty<Real>("gamma_damaged")),
    _eps_p(declareADProperty<RankTwoTensor>("eps_p")),
    _eps_e(declareADProperty<RankTwoTensor>("eps_e")),
    _sigma_d(declareADProperty<RankTwoTensor>("sigma_d")),
    _initial_damage(getADMaterialProperty<Real>("initial_damage"))
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

  //set initial damage
  _alpha_damagedvar[_qp] = _initial_damage[_qp];

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