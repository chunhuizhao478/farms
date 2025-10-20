//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Kernel.h"

/**
 * Diffusion operator for the breakage evolution equation.
 * Implements D * grad(B) · grad(phi_i) with an optional taper (1 - B)
 * handled directly in the residual to keep the front smooth as B -> 1.
 */
class BreakageEvolutionDiffusion : public Kernel
{
public:
  static InputParameters validParams();

  BreakageEvolutionDiffusion(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;

  virtual Real computeQpJacobian() override;

  const MaterialProperty<Real> & _D_diffusion;
};
