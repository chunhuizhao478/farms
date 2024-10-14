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
#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include "RotationTensor.h"
#include "DerivativeMaterialInterface.h"

/**
 * ComputeGeneralDamageBreakageStressBase3D is the direct base class for stress calculator materials that may
 * leverage quantities based on the displaced mesh (like the UMAT plugins) rather than solely using
 * strain tensors computed by separate MOOSE material objects (those classes should directly derive
 * from ComputeStressBase, which in turn directly derives from ComputeGeneralDamageBreakageStressBase3D).
 */
class ComputeGeneralDamageBreakageStressBase3D : public DerivativeMaterialInterface<Material>
{
public:
  static InputParameters validParams();

  ComputeGeneralDamageBreakageStressBase3D(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  /**
   * Compute the stress and store it in the _stress material property
   * for the current quadrature point
   **/
  virtual void computeQpStress() = 0;

  /// Function: Compute initial strain based on initial stress
  /// void computeInitialStrain();

  /// Base name prepended to all material property names to allow for
  /// multi-material systems
  const std::string _base_name;

  /// Mechanical strain material property
  const MaterialProperty<RankTwoTensor> & _mechanical_strain;
  /// Stress material property
  MaterialProperty<RankTwoTensor> & _stress;
  /// Elastic strain material property
  MaterialProperty<RankTwoTensor> & _elastic_strain;

  /// Extra stress tensor
  // const MaterialProperty<RankTwoTensor> & _extra_stress;

  /// initial stress components
  std::vector<const Function *> _initial_stress_fcn;

  /// derivative of stress w.r.t. strain (_dstress_dstrain)
  MaterialProperty<RankFourTensor> & _Jacobian_mult;

  /// additional material property
  /// alpha (damage parameter)
  MaterialProperty<Real> & _alpha_damagedvar;
  /// B (breakage parameter)
  MaterialProperty<Real> & _B;
  /// xi (strain invariant ratio)
  MaterialProperty<Real> & _xi;
  /// I1 (first strain invariant)
  MaterialProperty<Real> & _I1;
  /// I2 (second strain invariant)
  MaterialProperty<Real> & _I2;
  /// lambda (first lame const)
  MaterialProperty<Real> & _lambda;
  /// mu (shear modulus)
  MaterialProperty<Real> & _shear_modulus;
  /// gamma_damaged (damage modulus)
  MaterialProperty<Real> & _gamma_damaged;
  /// viscoelastic strain
  MaterialProperty<RankTwoTensor> & _eps_p;
  /// elastic strain
  MaterialProperty<RankTwoTensor> & _eps_e;
  /// total strain
  MaterialProperty<RankTwoTensor> & _eps_total;
  /// total stress tensor
  MaterialProperty<RankTwoTensor> & _sts_total;
  /// deviatoric stress tensor
  MaterialProperty<RankTwoTensor> & _sigma_d;
  /// eq strain rate
  MaterialProperty<Real> & _epsilon_eq;
  /// take initial value 
  /// lambda (first lame const)
  Real _lambda_o;
  /// mu (shear modulus)
  Real _shear_modulus_o;
};
