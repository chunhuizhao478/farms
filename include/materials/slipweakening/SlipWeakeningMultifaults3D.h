/* Debug */
#pragma once

#include "CZMComputeLocalTractionTotalBaseLSW3D.h"

class SlipWeakeningMultifaults3D : public CZMComputeLocalTractionTotalBaseLSW3D
{
public:
  static InputParameters validParams();
  SlipWeakeningMultifaults3D(const InputParameters & parameters);

protected:
  /// method computing the total traction and its derivatives
  void computeInterfaceTractionAndDerivatives() override;

  Real _Dc;

  const VariableValue & _nodal_area;

  const MaterialProperty<Real> & _density;

  const MaterialProperty<RankTwoTensor> & _rot;

  const VariableValue & _disp_slipweakening_x;
  const VariableValue & _disp_slipweakening_neighbor_x;
  const VariableValue & _disp_slipweakening_y;
  const VariableValue & _disp_slipweakening_neighbor_y;
  const VariableValue & _disp_slipweakening_z;
  const VariableValue & _disp_slipweakening_neighbor_z;

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

  const VariableValue & _mu_s;
  const VariableValue & _mu_d;

  const MaterialProperty<RankTwoTensor> & _sts_init;

  const bool _Co_coupled;
  const bool _T_coupled;

  const VariableValue * const _Co;
  const VariableValue * const _T;

  const MaterialProperty<Real> & _accumulated_slip_along_normal_old;
  const MaterialProperty<Real> & _accumulated_slip_along_strike_old;
  const MaterialProperty<Real> & _accumulated_slip_along_dip_old;

  const MaterialProperty<Real> & _slip_along_normal_old;
  const MaterialProperty<Real> & _slip_along_strike_old;
  const MaterialProperty<Real> & _slip_along_dip_old;

};