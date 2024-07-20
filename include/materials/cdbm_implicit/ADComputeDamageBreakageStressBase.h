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
  ADMaterialProperty<Real> & _B_breakagevar;
  /// xi (strain invariant ratio)
  ADMaterialProperty<Real> & _xi;
  /// I1 (first strain invariant)
  ADMaterialProperty<Real> & _I1;
  /// I2 (second strain invariant)
  ADMaterialProperty<Real> & _I2;
  /// mu (shear modulus)
  ADMaterialProperty<Real> & _shear_modulus;
  /// gamma_damaged (damage modulus)
  ADMaterialProperty<Real> & _gamma_damaged;
  /// viscoelastic strain
  ADMaterialProperty<RankTwoTensor> & _eps_p;
  /// elastic strain
  ADMaterialProperty<RankTwoTensor> & _eps_e;
  /// deviatoric stress
  ADMaterialProperty<RankTwoTensor> & _sigma_d;

};