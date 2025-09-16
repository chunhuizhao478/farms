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

class SlipWeakeningFrictionczm3dCDBM : public CZMComputeLocalTractionTotalBase
{
public:
  static InputParameters validParams();
  SlipWeakeningFrictionczm3dCDBM(const InputParameters & parameters);

protected:
  /// method computing the total traction and its derivatives
  void computeInterfaceTractionAndDerivatives() override;
  
  Real _mu_s;
  Real _mu_d;
  Real _Dc;
  Real _len;

  const MaterialProperty<Real> & _density;

  const MaterialProperty<RankTwoTensor> & _rot;

  const VariableValue & _disp_slipweakening_x;
  const VariableValue & _disp_slipweakening_neighbor_x;
  const VariableValue & _disp_slipweakening_y;
  const VariableValue & _disp_slipweakening_neighbor_y;
  const VariableValue & _disp_slipweakening_z;
  const VariableValue & _disp_slipweakening_neighbor_z;

  const VariableValue & _vel_slipweakening_x;
  const VariableValue & _vel_slipweakening_neighbor_x;
  const VariableValue & _vel_slipweakening_y;
  const VariableValue & _vel_slipweakening_neighbor_y;
  const VariableValue & _vel_slipweakening_z;
  const VariableValue & _vel_slipweakening_neighbor_z;

  const VariableValue & _reaction_slipweakening_x;
  const VariableValue & _reaction_slipweakening_neighbor_x;
  const VariableValue & _reaction_slipweakening_y;
  const VariableValue & _reaction_slipweakening_neighbor_y;
  const VariableValue & _reaction_slipweakening_z;
  const VariableValue & _reaction_slipweakening_neighbor_z;

  const VariableValue & _disp_slipweakening_x_old;
  const VariableValue & _disp_slipweakening_neighbor_x_old;
  const VariableValue & _disp_slipweakening_y_old;
  const VariableValue & _disp_slipweakening_neighbor_y_old;
  const VariableValue & _disp_slipweakening_z_old;
  const VariableValue & _disp_slipweakening_neighbor_z_old;

  // Define material properties for slip/slip rate
  MaterialProperty<Real> & _displacement_jump_strike; 
  MaterialProperty<Real> & _displacement_jump_dip;
  MaterialProperty<Real> & _displacement_jump_normal;
  MaterialProperty<Real> & _displacement_jump_rate_strike;
  MaterialProperty<Real> & _displacement_jump_rate_dip;
  MaterialProperty<Real> & _displacement_jump_rate_normal;

  // Define material properties for collecting total shear/normal traction
  MaterialProperty<Real> & _traction_strike;
  MaterialProperty<Real> & _traction_dip;
  MaterialProperty<Real> & _traction_normal;

  // Initial shear stress tensor
  const MaterialProperty<RankTwoTensor> & _static_initial_stress_tensor;

  // Use forced rupture
  bool _use_forced_rupture;
  Real _t0;
  const VariableValue & _cohesion_aux;
  const VariableValue & _forced_rupture_aux;
  const VariableValue & _fluid_pressure_aux;
};