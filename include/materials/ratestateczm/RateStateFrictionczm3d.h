//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "CZMComputeLocalTractionTotalBase.h"

/**
 * Rate-and-state friction with fast-velocity-weakening law for cohesive zone models
 * Based on the formulation in Premus et al. (2020)
 */
class RateStateFrictionczm3d : public CZMComputeLocalTractionTotalBase
{
public:
  static InputParameters validParams();

  RateStateFrictionczm3d(const InputParameters & parameters);

protected:
  void computeInterfaceTractionAndDerivatives() override;

  /// Helper functions for rate-and-state friction
  Real computeFSS(Real slip_rate);
  Real computePsiSS(Real slip_rate);
  Real updateStateVariable(Real slip_rate_current, Real psi_old);
  Real computeFrictionStrength(Real slip_rate, Real psi, Real normal_stress);
  Real solveNewtonForSlipRate(Real s_tilde, Real C, Real psi);

  // Rate-and-state friction parameters
  const Real _a;         ///< Direct effect parameter
  const Real _b;         ///< Evolution effect parameter
  const Real _L;         ///< Characteristic length scale
  const Real _f0;        ///< Reference friction coefficient
  const Real _fw;        ///< Weakened friction coefficient
  const Real _fLV;       ///< Low velocity friction coefficient
  const Real _s0;        ///< Reference slip rate
  const Real _sw;        ///< Weakening slip rate
  const Real _s_ini;     ///< Initial slip rate for numerical regularization

  // Geometric parameters
  const Real _len;       ///< Element edge length

  // Material properties
  const MaterialProperty<Real> & _density;
  const MaterialProperty<RankTwoTensor> & _rot;

  // Displacement variables
  const VariableValue & _disp_slipweakening_x;
  const VariableValue & _disp_slipweakening_neighbor_x;
  const VariableValue & _disp_slipweakening_y;
  const VariableValue & _disp_slipweakening_neighbor_y;
  const VariableValue & _disp_slipweakening_z;
  const VariableValue & _disp_slipweakening_neighbor_z;

  // Velocity variables
  const VariableValue & _vel_slipweakening_x;
  const VariableValue & _vel_slipweakening_neighbor_x;
  const VariableValue & _vel_slipweakening_y;
  const VariableValue & _vel_slipweakening_neighbor_y;
  const VariableValue & _vel_slipweakening_z;
  const VariableValue & _vel_slipweakening_neighbor_z;

  // Reaction variables
  const VariableValue & _reaction_slipweakening_x;
  const VariableValue & _reaction_slipweakening_neighbor_x;
  const VariableValue & _reaction_slipweakening_y;
  const VariableValue & _reaction_slipweakening_neighbor_y;
  const VariableValue & _reaction_slipweakening_z;
  const VariableValue & _reaction_slipweakening_neighbor_z;

  // Old displacement values
  const VariableValue & _disp_slipweakening_x_old;
  const VariableValue & _disp_slipweakening_neighbor_x_old;
  const VariableValue & _disp_slipweakening_y_old;
  const VariableValue & _disp_slipweakening_neighbor_y_old;
  const VariableValue & _disp_slipweakening_z_old;
  const VariableValue & _disp_slipweakening_neighbor_z_old;

  // State variable
  const VariableValue & _state_variable;
  const VariableValue & _state_variable_old;

  // Output material properties
  MaterialProperty<Real> & _displacement_jump_strike;
  MaterialProperty<Real> & _displacement_jump_dip;
  MaterialProperty<Real> & _displacement_jump_normal;
  MaterialProperty<Real> & _displacement_jump_rate_strike;
  MaterialProperty<Real> & _displacement_jump_rate_dip;
  MaterialProperty<Real> & _displacement_jump_rate_normal;
  MaterialProperty<Real> & _traction_strike;
  MaterialProperty<Real> & _traction_dip;
  MaterialProperty<Real> & _traction_normal;
  MaterialProperty<Real> & _slip_rate_magnitude;
  MaterialProperty<Real> & _friction_coefficient;
  MaterialProperty<Real> & _state_variable_updated;

  // Initial stress and auxiliary variables
  const MaterialProperty<RankTwoTensor> & _static_initial_stress_tensor;
  const VariableValue & _cohesion_aux;
  const VariableValue & _fluid_pressure_aux;

  // Newton solver parameters
  const Real _newton_tolerance;
  const unsigned int _newton_max_iterations;
};