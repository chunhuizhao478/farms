//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

/*
Material Description of Slip Weakening Friction Law 2D
*/

#pragma once

#include "CZMComputeLocalTractionTotalBase.h"

class SlipWeakeningFrictionczm2dCDBM : public CZMComputeLocalTractionTotalBase
{
public:
  static InputParameters validParams();
  SlipWeakeningFrictionczm2dCDBM(const InputParameters & parameters);

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

  const VariableValue & _reaction_slipweakening_x;
  const VariableValue & _reaction_slipweakening_neighbor_x;
  const VariableValue & _reaction_slipweakening_y;
  const VariableValue & _reaction_slipweakening_neighbor_y;

  const VariableValue & _disp_slipweakening_x_old;
  const VariableValue & _disp_slipweakening_neighbor_x_old;
  const VariableValue & _disp_slipweakening_y_old;
  const VariableValue & _disp_slipweakening_neighbor_y_old;

  // Define material properties for slip/slip rate
  MaterialProperty<Real> & _displacement_jump_strike; 
  MaterialProperty<Real> & _displacement_jump_normal;
  MaterialProperty<Real> & _displacement_jump_rate_strike;
  MaterialProperty<Real> & _displacement_jump_rate_normal;

  // Define material properties for collecting total shear/normal traction
  MaterialProperty<Real> & _traction_strike;
  MaterialProperty<Real> & _traction_normal;

  // Initial shear stress tensor
  const MaterialProperty<RankTwoTensor> & _static_initial_stress_tensor;

  /// peak shear stress
  Real _peak_shear_stress;
  /// nucleation center
  std::vector<Real> _nucl_center;
  /// nucleation radius
  Real _nucl_radius;

};
