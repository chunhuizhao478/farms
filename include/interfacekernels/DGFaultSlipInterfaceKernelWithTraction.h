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
 * DGFaultSlipInterfaceKernelWithTraction enforces prescribed slip at a fault interface
 * using DG/Nitsche method AND outputs the computed traction to an AuxVariable.
 *
 * This ensures the traction is computed using EXACTLY the same formulas and values
 * as the slip BC enforcement, avoiding any inconsistency.
 *
 * The traction formula follows Tandem's convention:
 *   τ = μ * {{∂u/∂n}} - κ * ([[u]] - slip_prescribed) + τ₀
 *
 * where:
 *   - {{∂u/∂n}} = 0.5 * (∇u_elem + ∇u_neighbor) · n
 *   - [[u]] = u_neighbor - u_elem (MOOSE convention)
 *   - κ = max(σ*μ/h, penalty/h)
 *   - τ₀ = pre-stress
 */
class DGFaultSlipInterfaceKernelWithTraction : public InterfaceKernel
{
public:
  static InputParameters validParams();

  DGFaultSlipInterfaceKernelWithTraction(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual(Moose::DGResidualType type) override;
  virtual Real computeQpJacobian(Moose::DGJacobianType type) override;

  /// Compute traction at current qp (called during residual computation)
  Real computeTraction() const;

  /// Prescribed slip (from SEASSlipAux)
  const VariableValue & _slip_prescribed;

  /// Shear modulus for consistent flux
  const MaterialProperty<Real> & _shear_modulus;

  /// Penalty parameter
  const Real _penalty;

  /// SIPG parameter (epsilon)
  const Real _epsilon;

  /// Symmetry parameter (sigma)
  const Real _sigma;

  /// Pre-stress to add to traction
  const Real _tau_pre;

  /// Writable traction AuxVariable
  MooseVariable * _traction_var;
};
