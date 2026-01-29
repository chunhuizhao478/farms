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
 * DGElasticTractionMaterial computes the elastic shear traction on a fault
 * interface for SEAS (Sequences of Earthquakes and Aseismic Slip) simulations.
 *
 * Three modes are available:
 *
 * 1. "gradient" mode (default): Uses simple DG formulation
 *    τ = μ * {{∂u/∂n}} = μ * (0.5 * (∇u_elem + ∇u_neighbor)) · n
 *    Note: This mode may produce incorrect stress rates for quasi-static SEAS
 *    problems because it ignores the penalty correction term.
 *
 * 2. "stiffness" mode: Uses discrete stiffness relationship
 *    τ = K * (Vp*t - slip)
 *    where K = μ / fault_depth is the stiffness for a fault in half-space.
 *    This mode correctly captures the quasi-static stress evolution.
 *
 * 3. "dg_consistent" mode: Uses full DG-consistent formula (Tandem-style)
 *    τ = μ * {{∂u/∂n}} - κ * ([[u]] - slip_prescribed)
 *    where κ = max(σ*μ/h, penalty/h) is the penalty coefficient.
 *    This is the correct DG formulation that accounts for the penalty
 *    enforcement of the slip boundary condition.
 */
class DGElasticTractionMaterial : public InterfaceMaterial
{
public:
  static InputParameters validParams();

  DGElasticTractionMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  /// Displacement on element side
  const VariableValue & _u;
  /// Displacement on neighbor side
  const VariableValue & _u_neighbor;

  /// Displacement gradient on element side
  const VariableGradient & _grad_u;
  /// Displacement gradient on neighbor side
  const VariableGradient & _grad_u_neighbor;

  /// Prescribed slip variable
  const VariableValue & _slip_prescribed;

  /// Shear modulus
  const Real _shear_modulus;

  /// DG penalty scaling factor
  const Real _sigma;

  /// Explicit penalty parameter
  const Real _penalty;

  /// Traction scaling factor
  const Real _traction_scale;

  /// Current element volume
  const Real & _current_elem_volume;

  /// Current side/face volume
  const Real & _current_side_volume;

  /// Output traction property
  MaterialProperty<Real> & _elastic_traction;

  /// Mode for traction computation: "gradient" or "stiffness"
  const MooseEnum _traction_mode;

  /// Plate rate for stiffness mode (m/s)
  const Real _plate_rate;

  /// Fault depth for stiffness mode (m)
  const Real _fault_depth;

  /// Domain half-width for stiffness mode (m)
  const Real _domain_half_width;
};
