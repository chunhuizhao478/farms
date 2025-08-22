//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ElementIntegralPostprocessor.h"
#include "MooseVariableInterface.h"

/**
 * This postprocessor computes a volume integral of the specified variable.
 *
 * Note that specializations of this integral are possible by deriving from this
 * class and overriding computeQpIntegral().
 */
class FirstStepElementIntegralVariablePostprocessor : public ElementIntegralPostprocessor,
                                             public MooseVariableInterface<Real>
{
public:
  static InputParameters validParams();

  FirstStepElementIntegralVariablePostprocessor(const InputParameters & parameters);

protected:
  // ElementIntegralPostprocessor API
  Real computeQpIntegral() override;
  Real getValue() const override;

  // Freeze the value after first step by skipping base accumulation after capture
  void initialize() override;
  void execute() override;
  void finalize() override;

  /// Holds the solution at current quadrature points
  const VariableValue & _u;
  /// Holds the solution at previous quadrature points
  const VariableValue & _u_old;
  /// Holds the solution gradient at the current quadrature points
  const VariableGradient & _grad_u;
  /// Option to use absolute variable value
  bool _use_abs_value;

  // Stored value of the first-step integral and a flag indicating capture
  Real _first_value;
  bool _captured;
};
