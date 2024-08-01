//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "FarmsSlipWeakeningBase.h"

/**
 * Implement slip weakening friction law in CZM
 **/
class FarmsSlipWeakeningCZM : public FarmsSlipWeakeningBase
{
public:
  static InputParameters validParams();
  FarmsSlipWeakeningCZM(const InputParameters & parameters);
protected:
  virtual Real computeTractionAndDisplacements() override;
  virtual RealTensorValue computeTractionDerivatives() override;
  virtual void initQpStatefulProperties() override;
  
  Real _Dc;

  const MaterialProperty<Real> & _density;

  const MaterialProperty<RealTensorValue> & _rot;

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

  const VariableValue & _reaction_damp_x;
  const VariableValue & _reaction_damp_neighbor_x;
  const VariableValue & _reaction_damp_y;
  const VariableValue & _reaction_damp_neighbor_y;
  const VariableValue & _reaction_damp_z;
  const VariableValue & _reaction_damp_neighbor_z;

  const VariableValue & _disp_slipweakening_x_old;
  const VariableValue & _disp_slipweakening_neighbor_x_old;
  const VariableValue & _disp_slipweakening_y_old;
  const VariableValue & _disp_slipweakening_neighbor_y_old;
  const VariableValue & _disp_slipweakening_z_old;
  const VariableValue & _disp_slipweakening_neighbor_z_old;

  const VariableValue & _elem_length;

  const VariableValue & _mu_s;
  const VariableValue & _mu_d;

  const MaterialProperty<RankTwoTensor> & _sts_init;

  const bool _Co_coupled;
  const bool _T_coupled;

  const VariableValue * const _Co;
  const VariableValue * const _T;

  const MaterialProperty<RealVectorValue> & _displacements_plus_old;
  const MaterialProperty<RealVectorValue> & _displacements_minus_old;
  const MaterialProperty<RealVectorValue> & _displacements_plus_older;
  const MaterialProperty<RealVectorValue> & _displacements_minus_older;
  const MaterialProperty<RealVectorValue> & _velocities_plus_old;
  const MaterialProperty<RealVectorValue> & _velocities_minus_old;

  const MaterialProperty<Real> & _absolute_slip_old;

  const MaterialProperty<RealVectorValue> & _below_strength_marker_old; //this is used to retrieve the marker for all elements the traction is below strength

  const MaterialProperty<RealVectorValue> & _R_plus_local_vec_old;
  const MaterialProperty<RealVectorValue> & _R_minus_local_vec_old;

  const MaterialProperty<RealVectorValue> & _traction_total_local_old;
};
