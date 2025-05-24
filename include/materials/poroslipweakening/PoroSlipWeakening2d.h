/* 
Material Description of Slip Weakening Friction v3
Generalize the computation of sticking traction using consistent displacement jump and nodal reaction forces
*/

#pragma once

#include "PoroCZMComputeLocalTractionBase.h"

class PoroSlipWeakening2d : public PoroCZMComputeLocalTractionBase
{
public:
  static InputParameters validParams();
  PoroSlipWeakening2d(const InputParameters & parameters);

protected:
  /// method computing the total traction and its derivatives
  void computeInterfaceTractionAndDerivatives() override;

  Real _T2_o;
  Real _mu_d;
  Real _Dc;
  Real _nodal_area2;
  

  const MaterialProperty<Real> &_biot_coefficient;
  const MaterialProperty<RankTwoTensor> & _rot;
  const VariableValue & _nodal_area;
  const MaterialProperty<Real> & _rhof;
  const MaterialProperty<Real> & _density;

  const MaterialProperty<RankTwoTensor> & _stress;
  const MaterialProperty<RankTwoTensor> & _dstress;

  const std::string _permeability_type;

  const VariableValue & _interface_pressure_plus; 
  const VariableValue & _interface_pressure_minus;
  const VariableValue & _interface_pressure_plus_older; 
  const VariableValue & _interface_pressure_minus_older;

  const MaterialProperty<RealVectorValue> & _interface_displacement_jump_old;
  const MaterialProperty<RealVectorValue> & _interface_displacement_jump_older;

  const VariableValue & _reaction_x;
  const VariableValue & _reaction_neighbor_x;
  const VariableValue & _reaction_y;
  const VariableValue & _reaction_neighbor_y;

  const VariableValue & _reaction_damp_x;
  const VariableValue & _reaction_neighbor_damp_x;
  const VariableValue & _reaction_damp_y;
  const VariableValue & _reaction_neighbor_damp_y;

  const VariableValue & _reaction_pressure_x;
  const VariableValue & _reaction_neighbor_pressure_x;
  const VariableValue & _reaction_pressure_y;
  const VariableValue & _reaction_neighbor_pressure_y;

  const VariableValue & _jacobian_x;
  const VariableValue & _jacobian_neighbor_x;
  const VariableValue & _jacobian_y;
  const VariableValue & _jacobian_neighbor_y;

  const VariableValue & _jacobian_damp_x;
  const VariableValue & _jacobian_neighbor_damp_x;
  const VariableValue & _jacobian_damp_y;
  const VariableValue & _jacobian_neighbor_damp_y;

  const VariableValue & _jacobian_pressure_x;
  const VariableValue & _jacobian_neighbor_pressure_x;
  const VariableValue & _jacobian_pressure_y;
  const VariableValue & _jacobian_neighbor_pressure_y;

  const MaterialProperty<RealVectorValue> & _interface_fluid_vel_jump_old;

  const VariableValue & _fluid_disp_x;
  const VariableValue & _fluid_disp_neighbor_x;
  const VariableValue & _fluid_disp_y;
  const VariableValue & _fluid_disp_neighbor_y;


};