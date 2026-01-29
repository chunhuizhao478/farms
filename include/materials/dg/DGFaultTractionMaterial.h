//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "InterfaceMaterial.h"

/**
 * DGFaultTractionMaterial computes fault traction for DG SEAS simulations.
 * This is a simple material that provides either:
 * 1. A prescribed constant traction
 * 2. A linear slip-weakening or slip-strengthening traction
 *
 * For SEAS simulations, this would be replaced by a full rate-state friction material.
 */
class DGFaultTractionMaterial : public InterfaceMaterial
{
public:
  static InputParameters validParams();

  DGFaultTractionMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  /// Displacement variable on element side
  const VariableValue & _u;

  /// Displacement variable on neighbor side
  const VariableValue & _u_neighbor;

  /// Prescribed background traction
  const Real _tau0;

  /// Slip-traction coupling coefficient (dτ/ds)
  const Real _dtau_ds;

  /// Output: fault traction
  MaterialProperty<Real> & _fault_traction;

  /// Output: derivative of traction w.r.t. slip
  MaterialProperty<Real> & _dtraction_dslip;
};
