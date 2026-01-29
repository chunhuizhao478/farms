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
 * DGFaultSlipInterfaceKernel enforces a prescribed slip (displacement jump)
 * at a fault interface using a penalty method.
 *
 * This is part of the staggered SEAS solver where:
 * - Slip is computed explicitly from friction (AuxKernels)
 * - Elasticity is solved implicitly with prescribed slip
 *
 * The constraint enforced is:
 *   [u] = u_elem - u_neighbor = slip_prescribed
 *
 * Using Nitsche's method / penalty:
 *   penalty * ([u] - slip_prescribed) = 0
 *
 * This makes the elasticity problem LINEAR because slip is prescribed.
 */
class DGFaultSlipInterfaceKernel : public InterfaceKernel
{
public:
  static InputParameters validParams();

  DGFaultSlipInterfaceKernel(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual(Moose::DGResidualType type) override;
  virtual Real computeQpJacobian(Moose::DGJacobianType type) override;

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
};
