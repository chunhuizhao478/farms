//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "DGKernel.h"

/**
 * DGElasticityAntiplane implements the Symmetric Interior Penalty Galerkin (SIPG)
 * method for 2D antiplane shear elasticity. This is based on the formulation
 * in the Tandem paper for SEAS simulations.
 *
 * The strong form: -mu * Laplacian(w) = 0
 *
 * The SIPG bilinear form for internal facets:
 * a(w,v) = -{{mu*grad(w)}}*[[v]] - epsilon*{{mu*grad(v)}}*[[w]] + sigma/h*[[w]]*[[v]]
 *
 * where:
 *   {{.}} = average operator
 *   [[.]] = jump operator
 *   epsilon = 1 (SIPG), -1 (NIPG), 0 (IIPG)
 *   sigma = penalty parameter
 */
class DGElasticityAntiplane : public DGKernel
{
public:
  static InputParameters validParams();

  DGElasticityAntiplane(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual(Moose::DGResidualType type) override;
  virtual Real computeQpJacobian(Moose::DGJacobianType type) override;

  /// Symmetry parameter: 1 = SIPG, -1 = NIPG, 0 = IIPG
  const Real _epsilon;

  /// Penalty parameter multiplier
  const Real _sigma;

  /// Shear modulus on element side
  const MaterialProperty<Real> & _mu;

  /// Shear modulus on neighbor side
  const MaterialProperty<Real> & _mu_neighbor;
};
