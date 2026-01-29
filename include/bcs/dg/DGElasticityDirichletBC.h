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
 * DGElasticityDirichletBC implements a Dirichlet boundary condition for
 * DG elasticity using the SIPG method. This enforces the displacement weakly
 * using interior penalty terms at the boundary.
 *
 * For a boundary with prescribed displacement g:
 *   R = -mu*grad(u).n * v + epsilon*(u-g)*mu*grad(v).n + sigma/h*(u-g)*v
 */
class DGElasticityDirichletBC : public IntegratedBC
{
public:
  static InputParameters validParams();

  DGElasticityDirichletBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  /// Prescribed displacement value
  const Real _value;

  /// Optional function for prescribed displacement
  const Function * const _func;

  /// Symmetry parameter: 1 = SIPG, -1 = NIPG, 0 = IIPG
  const Real _epsilon;

  /// Penalty parameter multiplier
  const Real _sigma;

  /// Shear modulus
  const MaterialProperty<Real> & _mu;

  /// Helper to get prescribed value
  Real prescribedValue() const;
};
