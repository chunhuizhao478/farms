//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ADComputeStressBase.h"

/**
 * ADComputeDamageStressStaticDistributionDynamicCDBM computes the stress following linear elasticity theory (small strains)
 */
class ADComputeDamageStressStaticDistributionDynamicCDBM : public ADComputeStressBase
{
public:
  static InputParameters validParams();

  ADComputeDamageStressStaticDistributionDynamicCDBM(const InputParameters & parameters);

  virtual void initialSetup() override;

protected:
  virtual void computeQpStress() override;

  /// @brief Compute gamma_r
  /// @return gamma_r
  ADReal computegammar();

  /// @brief Compute breakage coefficients
  /// @param gamma_damaged_r
  /// @return a0 a1 a2 a3 in a vector
  std::vector<ADReal> computecoefficients(ADReal gamma_damaged_r);  

  /// @brief Compute first root of hessian matrix
  /// @param xi 
  /// @return the first root of critical alpha_cr
  ADReal alphacr_root1(ADReal xi, ADReal gamma_damaged_r);  

  /// Material property initial damage profile

  /// initial lambda value 
  ADReal _lambda_o;

  /// initial shear modulus value
  ADReal _shear_modulus_o;

  /// xi_o value
  ADReal _xi_o;

  /// xi_d value
  ADReal _xi_d;

  /// chi value
  ADReal _chi;

  /// @brief initial damage value
  const ADMaterialProperty<Real> & _initial_damage_val;

  /// @brief initial breakage value
  const ADMaterialProperty<Real> & _initial_breakage_val;

};
