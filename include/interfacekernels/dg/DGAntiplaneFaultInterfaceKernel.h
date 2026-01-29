//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "InterfaceKernel.h"

/**
 * DGAntiplaneFaultInterfaceKernel implements the fault interface condition for
 * DG antiplane shear elasticity in SEAS simulations.
 *
 * Based on the Tandem paper numerical flux table for fault facets:
 *   û = u - s/2 (on + side), û = u + s/2 (on - side)
 *   σ̂·n = τ (traction from friction law)
 *
 * The weak form contribution from fault interface:
 *   ∫_Γf τ * [[v]] dA
 *
 * where τ is the fault traction (from rate-state friction or prescribed),
 * and [[v]] = v+ - v- is the jump in test function.
 *
 * The fault traction comes from a material property computed by the
 * rate-state friction material.
 */
class DGAntiplaneFaultInterfaceKernel : public InterfaceKernel
{
public:
  static InputParameters validParams();

  DGAntiplaneFaultInterfaceKernel(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual(Moose::DGResidualType type) override;
  virtual Real computeQpJacobian(Moose::DGJacobianType type) override;

  /// Fault traction from friction law (material property)
  const MaterialProperty<Real> & _fault_traction;

  /// Derivative of traction w.r.t. slip (for Jacobian)
  const MaterialProperty<Real> * _dtraction_dslip;

  /// Derivative of traction w.r.t. slip rate (for Jacobian)
  const MaterialProperty<Real> * _dtraction_dslip_rate;

  /// Whether to use penalty enforcement for slip constraint
  const bool _use_penalty;

  /// Penalty parameter for slip enforcement (if used)
  const Real _penalty;

  /// Shear modulus for penalty calculation
  const MaterialProperty<Real> * _mu;
  const MaterialProperty<Real> * _mu_neighbor;
};
