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
 * This kernel implements the Laplacian operator for damage evolution diffusion equation:
 * $(1-B) D \nabla \alpha \cdot \nabla \phi_i$
 */
class DamageEvolutionDiffusion : public Kernel
{
public:
  static InputParameters validParams();

  DamageEvolutionDiffusion(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;

  virtual Real computeQpJacobian() override;

  const MaterialProperty<Real> & _B_breakagevar; //breakage variable
  const MaterialProperty<Real> & _D_diffusion; //diffusion coefficient
};