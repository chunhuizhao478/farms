//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

/*
Material Description of Slip Weakening Friction 3d
*/

#pragma once

#include "CZMComputeLocalTractionTotalBase.h"

class RateStateFrictionczm3d : public CZMComputeLocalTractionTotalBase
{
public:
  static InputParameters validParams();
  RateStateFrictionczm3d(const InputParameters & parameters);

protected:
  /// method computing the total traction and its derivatives
  void computeInterfaceTractionAndDerivatives() override;

  void initQpStatefulProperties() override;

  /// declare material properties
  /// state variable (scalar)
  MaterialProperty<Real> & _statevar;
  const MaterialProperty<Real> & _statevar_old;

  /// slip rate (vector)
  MaterialProperty<RealVectorValue> & _slipratevar;
  const MaterialProperty<RealVectorValue> & _slipratevar_old;

  Real _T1_o;
  Real _T2_o;
  Real _T3_o;

  //rate-and-state friction coefficients
  Real _len;
  Real _f_o;
  Real _rsf_a;        // constant value (used if rsf_a_var not provided)
  Real _rsf_b;
  Real _rsf_L;
  Real _delta_o;
  Real _statevar_init; // constant value (used if statevar_init_var not provided)
  Real _sliprate_strike_init;

  // Spatially variable RSF parameters (optional coupled variables)
  const bool _use_coupled_rsf_a;
  const VariableValue * _rsf_a_var;

  const bool _use_coupled_statevar_init;
  const VariableValue * _statevar_init_var;

  const MaterialProperty<Real> & _density;

  const MaterialProperty<RankTwoTensor> & _rot;

  const VariableValue & _disp_x;
  const VariableValue & _disp_neighbor_x;
  const VariableValue & _disp_y;
  const VariableValue & _disp_neighbor_y;
  const VariableValue & _disp_z;
  const VariableValue & _disp_neighbor_z;

  const VariableValue & _vel_x;
  const VariableValue & _vel_neighbor_x;
  const VariableValue & _vel_y;
  const VariableValue & _vel_neighbor_y;
  const VariableValue & _vel_z;
  const VariableValue & _vel_neighbor_z;

  const VariableValue & _reaction_x;
  const VariableValue & _reaction_neighbor_x;
  const VariableValue & _reaction_y;
  const VariableValue & _reaction_neighbor_y;
  const VariableValue & _reaction_z;
  const VariableValue & _reaction_neighbor_z;

  const VariableValue & _disp_x_old;
  const VariableValue & _disp_neighbor_x_old;
  const VariableValue & _disp_y_old;
  const VariableValue & _disp_neighbor_y_old;
  const VariableValue & _disp_z_old;
  const VariableValue & _disp_neighbor_z_old;

  //shear stress perturbation
  ///Measure from current time step
  const VariableValue & _Ts_perturb;
  const VariableValue & _Ts_perturb_old;
};
