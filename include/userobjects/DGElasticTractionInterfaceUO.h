//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "InterfaceQpUserObjectBase.h"

/**
 * DGElasticTractionInterfaceUO computes the DG-consistent elastic traction
 * at fault interfaces for SEAS problems using the SIPG formulation.
 *
 * Following Tandem's Poisson operator:
 *   traction = mu * {{grad_u}} . n - penalty * slip_error
 *
 * where:
 *   {{grad_u}} = 0.5 * (grad_u_elem + grad_u_neighbor) is the average gradient
 *   n = unit normal from element to neighbor
 *   slip_error = (u_neighbor - u_elem) - slip_prescribed
 *   penalty = max(sigma * mu / h, explicit_penalty / h)
 *
 * This UserObject is designed for the MultiApp SEAS architecture where
 * traction needs to be transferred between apps.
 */
class DGElasticTractionInterfaceUO : public InterfaceQpUserObjectBase
{
public:
  static InputParameters validParams();

  DGElasticTractionInterfaceUO(const InputParameters & parameters);

protected:
  virtual Real computeRealValue(const unsigned int qp) override;

  /// Displacement on element side
  const VariableValue & _u;
  /// Displacement on neighbor side
  const VariableValue & _u_neighbor;

  /// Displacement gradient on element side
  const VariableGradient & _grad_u;
  /// Displacement gradient on neighbor side
  const VariableGradient & _grad_u_neighbor;

  /// Prescribed slip variable (optional)
  const VariableValue * _slip_prescribed;

  /// Shear modulus
  const Real _shear_modulus;

  /// DG penalty scaling factor
  const Real _sigma;

  /// Explicit penalty parameter
  const Real _penalty;

  /// Pre-stress/initial traction
  const Real _tau_pre;

  /// Reference to assembly for volume computation
  Assembly & _assembly;
};
