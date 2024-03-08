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

  //compute initial strain based on initial stress
  /// lambda (first lame const)
  _lambda[_qp] = _lambda_o;
  /// mu (shear modulus)
  _shear_modulus[_qp] = _shear_modulus_o;
  /// gamma_damaged (damage modulus)
  _gamma_damaged[_qp] = 0.0;

  //allpha, B
  _alpha_damagedvar[_qp] = 0.0;
  _B[_qp] = 0.0;

  //Convert (lambda_o,shear_modulus_o) to (youngs_modulus_o,poisson_ratio_o)
  ADReal youngs_modulus_o = _shear_modulus_o * ( 3 * _lambda_o + 2 * _shear_modulus_o ) / ( _lambda_o + _shear_modulus_o );
  ADReal poisson_ratio_o = _lambda_o / ( 2 * ( _lambda_o + _shear_modulus_o ));

  //Get stress components
  ADRankTwoTensor stress_initial = _static_initial_stress_tensor[_qp];

  //Get stress components
  ADReal sts11_init = stress_initial(0,0);
  ADReal sts12_init = stress_initial(0,1);
  ADReal sts22_init = stress_initial(1,1);
  //ADReal sts33_init = stress_initial(2,2);
  ADReal sts33_init = poisson_ratio_o * ( sts11_init + sts22_init );

  //Compute strain components using Hooke's Law
  ADReal eps11_init = 1.0 / youngs_modulus_o * ( sts11_init - poisson_ratio_o * ( sts22_init + sts33_init ) );
  ADReal eps22_init = 1.0 / youngs_modulus_o * ( sts22_init - poisson_ratio_o * ( sts11_init + sts33_init ) ); 
  ADReal eps12_init = 1.0 / youngs_modulus_o * ( ( 1 + poisson_ratio_o ) * sts12_init                       );
  ADReal eps13_init = 0.0;
  ADReal eps23_init = 0.0;
  ADReal eps33_init = 0.0;

  //Compute xi, I1, I2
  ADReal I1_init = eps11_init + eps22_init + eps33_init;
  ADReal I2_init = eps11_init * eps11_init + eps22_init * eps22_init + eps33_init * eps33_init + 2 * eps12_init * eps12_init + 2 * eps13_init * eps13_init + 2 * eps23_init * eps23_init;
  ADReal xi_init = I1_init / std::sqrt( I2_init );

  //Compute eps
  //eps_p
  _eps_p[_qp](0,0) = 0.0; _eps_p[_qp](0,1) = 0.0; _eps_p[_qp](0,2) = 0.0;
  _eps_p[_qp](1,0) = 0.0; _eps_p[_qp](1,1) = 0.0; _eps_p[_qp](1,2) = 0.0;
  _eps_p[_qp](2,0) = 0.0; _eps_p[_qp](2,1) = 0.0; _eps_p[_qp](2,2) = 0.0;
  //eps_e
  _eps_e[_qp](0,0) = eps11_init; _eps_e[_qp](0,1) = eps12_init; 
  _eps_e[_qp](1,0) = eps12_init; _eps_e[_qp](1,1) = eps22_init;
  _eps_e[_qp](2,2) = 0.0;
  //eps_total
  _eps_total[_qp](0,0) = eps11_init; _eps_total[_qp](0,1) = eps12_init;
  _eps_total[_qp](1,0) = eps12_init; _eps_total[_qp](1,1) = eps22_init;
  _eps_total[_qp](2,2) = 0.0;
  //sts_total
  _sts_total[_qp] = stress_initial;

  //I1
  _I1[_qp] = I1_init;
  //I2
  _I2[_qp] = I2_init;
  //xi
  _xi[_qp] = xi_init;  

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