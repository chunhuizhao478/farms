//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Kernel.h"

/**
 * StateEvolutionKernel implements the state evolution equation for rate-and-state friction.
 *
 * Aging Law: dθ/dt = 1 - V*θ/Dc
 *
 * This kernel provides the RHS contribution: -(1 - V*θ/Dc)
 * When combined with TimeDerivative kernel, gives:
 *   dθ/dt - (1 - V*θ/Dc) = 0  =>  dθ/dt = 1 - V*θ/Dc
 *
 * This is designed for use in the Friction SubApp of the MultiApp SEAS architecture,
 * where θ is solved implicitly while V is provided as an AuxVariable.
 *
 * Alternative: Slip Law: dθ/dt = -V*θ/Dc * ln(V*θ/Dc)
 * Can be selected via the 'evolution_law' parameter.
 */
class StateEvolutionKernel : public Kernel
{
public:
  static InputParameters validParams();

  StateEvolutionKernel(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  /// Slip rate V (from AuxVariable, updated by SlipRateAux)
  const VariableValue & _slip_rate;

  /// Coupled variable number for slip_rate (for off-diagonal Jacobian)
  const unsigned int _slip_rate_var;

  /// Critical slip distance Dc (m)
  const Real _Dc;

  /// Evolution law type: "aging" or "slip"
  const MooseEnum _evolution_law;
};
