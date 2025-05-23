/* 
Material Description of Slip Weakening Friction v3
Generalize the computation of sticking traction using consistent displacement jump and nodal reaction forces
*/

#pragma once

#include "CZMComputeLocalTractionTotalBase.h"

class PoroSlipWeakening2d_sub : public CZMComputeLocalTractionTotalBase
{
public:
  static InputParameters validParams();
  PoroSlipWeakening2d_sub(const InputParameters & parameters);

protected:
  /// method computing the total traction and its derivatives
  void computeInterfaceTractionAndDerivatives() override;

  Real _T2_o;
  Real _mu_d;
  Real _Dc;
  Real _nodal_area2;

  const MaterialProperty<RankTwoTensor> & _rot;
  const VariableValue & _nodal_area;
  const MaterialProperty<Real> & _rhof;
  const MaterialProperty<Real> & _density;

  const MaterialProperty<RankTwoTensor> & _stress;

  const std::string _permeability_type;

  const VariableValue & _interface_pressure_plus; 
  const VariableValue & _interface_pressure_minus;

  const VariableValue & _disp_x;
  const VariableValue & _disp_neighbor_x;
  const VariableValue & _disp_y;
  const VariableValue & _disp_neighbor_y;

  const VariableValue & _disp_x_older;
  const VariableValue & _disp_neighbor_x_older;
  const VariableValue & _disp_y_older;
  const VariableValue & _disp_neighbor_y_older;
  
  const VariableValue & _reaction_x;
  const VariableValue & _reaction_neighbor_x;
  const VariableValue & _reaction_y;
  const VariableValue & _reaction_neighbor_y;

  const VariableValue & _reaction_pressure_x;
  const VariableValue & _reaction_neighbor_pressure_x;
  const VariableValue & _reaction_pressure_y;
  const VariableValue & _reaction_neighbor_pressure_y;

  const VariableValue & _reaction_damp_x;
  const VariableValue & _reaction_neighbor_damp_x;
  const VariableValue & _reaction_damp_y;
  const VariableValue & _reaction_neighbor_damp_y;

  const VariableValue & _fluid_vel_x;
  const VariableValue & _fluid_vel_neighbor_x;
  const VariableValue & _fluid_vel_y;
  const VariableValue & _fluid_vel_neighbor_y;

  const VariableValue & _fluid_disp_x;
  const VariableValue & _fluid_disp_neighbor_x;
  const VariableValue & _fluid_disp_y;
  const VariableValue & _fluid_disp_neighbor_y;


};