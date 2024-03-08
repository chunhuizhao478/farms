//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Material.h"
#include "Function.h"
#include "ADRankTwoTensorForward.h"
#include "ADSymmetricRankTwoTensorForward.h"

// #define usingComputeDamageBreakageStressBaseMembers                                                              \
//   usingMaterialMembers;                                                                            \
//   using ADComputeDamageBreakageStressBase<RankTwoTensor>::_base_name;                                                  \
//   using ADComputeDamageBreakageStressBase<RankTwoTensor>::_mechanical_strain;                                          \
//   using ADComputeDamageBreakageStressBase<RankTwoTensor>::_stress;                                                     \
//   using ADComputeDamageBreakageStressBase<RankTwoTensor>::_elastic_strain;                                             \
//   using ADComputeDamageBreakageStressBase<RankTwoTensor>::_extra_stresses;                                             \
//   using ADComputeDamageBreakageStressBase<RankTwoTensor>::_initial_stress_fcn

/**
 * ADComputeDamageBreakageStressBase is the base class for stress tensors
 */
//template <typename RankTwoTensor>
class ADComputeDamageBreakageStressBase : public Material
{
public:
  static InputParameters validParams();

  ADComputeDamageBreakageStressBase(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties() override;
  virtual void computeQpProperties() override;
  virtual void computeQpStress() = 0;

  /// Base name of the material system
  const std::string _base_name;

  const ADMaterialProperty<RankTwoTensor> & _mechanical_strain;

  /// The stress tensor to be calculated
  ADMaterialProperty<RankTwoTensor> & _stress;
  ADMaterialProperty<RankTwoTensor> & _elastic_strain;

  /// Extra stress tensors
  std::vector<const MaterialProperty<RankTwoTensor> *> _extra_stresses;

  /// initial stress components
  std::vector<const Function *> _initial_stress_fcn;

  /// additional material property
  /// alpha (damage parameter)
  ADMaterialProperty<Real> & _alpha_damagedvar;
  /// B (breakage parameter)
  ADMaterialProperty<Real> & _B;
  /// xi (strain invariant ratio)
  ADMaterialProperty<Real> & _xi;
  /// I1 (first strain invariant)
  ADMaterialProperty<Real> & _I1;
  /// I2 (second strain invariant)
  ADMaterialProperty<Real> & _I2;
  /// lambda (first lame const)
  ADMaterialProperty<Real> & _lambda;
  /// mu (shear modulus)
  ADMaterialProperty<Real> & _shear_modulus;
  /// gamma_damaged (damage modulus)
  ADMaterialProperty<Real> & _gamma_damaged;
  /// viscoelastic strain
  ADMaterialProperty<RankTwoTensor> & _eps_p;
  /// elastic strain
  ADMaterialProperty<RankTwoTensor> & _eps_e;
  /// total strain
  ADMaterialProperty<RankTwoTensor> & _eps_total;
  /// total inital strain 
  ADMaterialProperty<RankTwoTensor> & _eps_total_init;
  /// total stress tensor
  ADMaterialProperty<RankTwoTensor> & _sts_total;
  /// initial stress tenosr
  const ADMaterialProperty<RankTwoTensor> & _static_initial_stress_tensor; 
  /// take initial value 
  /// lambda (first lame const)
  ADReal _lambda_o;
  /// mu (shear modulus)
  ADReal _shear_modulus_o; 

};

//typedef ADComputeDamageBreakageStressBase<RankTwoTensor> ADComputeDamageBreakageStressBase;