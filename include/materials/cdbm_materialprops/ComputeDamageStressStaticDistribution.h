//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ComputeStressBase.h"

/**
 * ComputeDamageStressStaticDistribution computes the stress following linear elasticity theory (small strains)
 */
class ComputeDamageStressStaticDistribution : public ComputeStressBase
{
public:
  static InputParameters validParams();

  ComputeDamageStressStaticDistribution(const InputParameters & parameters);

  virtual void initialSetup() override;

protected:
  virtual void computeQpStress() override;

  /// Material property initial damage profile

  /// initial lambda value 
  Real _lambda_o;

  /// initial shear modulus value
  Real _shear_modulus_o;

  /// xi_o value
  Real _xi_o;

  /// gamma_damage_r value
  Real _gamma_damage_r;

  const MaterialProperty<Real> & _initial_damage_val;

};
