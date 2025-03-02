/* 
This file aims to take care of special cases: fault opening/rejoining, reversal in the direction of slip in slip weakening simulation
Created by Chunhui Zhao, Dec 28, 2023
*/

#pragma once

#include "CZMComputeLocalTractionTotalBaseSWF2D.h"

class SlipWeakeningMultifaultsGeneralizedPropCSV : public CZMComputeLocalTractionTotalBaseSWF2D
{
public:
  static InputParameters validParams();
  SlipWeakeningMultifaultsGeneralizedPropCSV(const InputParameters & parameters);

protected:
  /// method computing the total traction and its derivatives
  void computeInterfaceTractionAndDerivatives() override;

  const VariableValue & _Dc;

  const VariableValue & _nodal_area;

  const MaterialProperty<Real> & _density;

  const MaterialProperty<RankTwoTensor> & _rot;

  const VariableValue & _disp_slipweakening_x;
  const VariableValue & _disp_slipweakening_neighbor_x;
  const VariableValue & _disp_slipweakening_y;
  const VariableValue & _disp_slipweakening_neighbor_y;

  const VariableValue & _reaction_slipweakening_x;
  const VariableValue & _reaction_slipweakening_neighbor_x;
  const VariableValue & _reaction_slipweakening_y;
  const VariableValue & _reaction_slipweakening_neighbor_y;

  const VariableValue & _disp_slipweakening_x_old;
  const VariableValue & _disp_slipweakening_neighbor_x_old;
  const VariableValue & _disp_slipweakening_y_old;
  const VariableValue & _disp_slipweakening_neighbor_y_old;

  const VariableValue & _mu_s;
  const VariableValue & _mu_d;
 
  const VariableValue & _ini_shear_sts;
  const VariableValue & _ini_normal_sts;

  const MaterialProperty<Real> & _flag_track_opening_old;
  const MaterialProperty<Real> & _flag_track_activecase_old;

  const MaterialProperty<Real> & _jump_track_opening_old;
  const MaterialProperty<Real> & _jump_track_reversal_old;

  const MaterialProperty<Real> & _T1_old;

};