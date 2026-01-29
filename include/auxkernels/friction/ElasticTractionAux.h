//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "AuxKernel.h"

/**
 * ElasticTractionAux computes the elastic shear traction from the displacement
 * gradient at a boundary.
 *
 * For antiplane shear (mode III):
 *   τ = μ * ∂w/∂n
 *
 * where n is the normal direction to the fault.
 *
 * This is designed for the Elasticity MainApp in the MultiApp SEAS architecture.
 * The computed traction is then transferred to the Friction SubApp.
 *
 * For general 3D elasticity, the traction would be:
 *   t = σ · n = C : ε(u) · n
 */
class ElasticTractionAux : public AuxKernel
{
public:
  static InputParameters validParams();

  ElasticTractionAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  /// Displacement variable
  const VariableValue & _disp;

  /// Displacement gradient
  const VariableGradient & _grad_disp;

  /// Shear modulus
  const Real _shear_modulus;

  /// Normal direction component (0=x, 1=y, 2=z)
  const unsigned int _normal_component;
};
