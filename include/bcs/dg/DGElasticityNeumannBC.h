//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "IntegratedBC.h"

/**
 * DGElasticityNeumannBC implements a Neumann (traction) boundary condition for
 * DG elasticity. For a boundary with prescribed traction t:
 *
 * According to the Tandem paper numerical flux table for Neumann boundaries:
 *   û = u  (displacement is continuous)
 *   σ̂·n = t  (traction is prescribed)
 *
 * Residual contribution: -t * v
 *
 * For traction-free (free surface), set traction = 0.
 */
class DGElasticityNeumannBC : public IntegratedBC
{
public:
  static InputParameters validParams();

  DGElasticityNeumannBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  /// Prescribed traction value
  const Real _traction;

  /// Optional function for prescribed traction
  const Function * const _func;
};
