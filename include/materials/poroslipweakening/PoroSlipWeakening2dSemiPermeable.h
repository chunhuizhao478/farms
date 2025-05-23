/* 
Material Description of Slip Weakening Friction v3
Generalize the computation of sticking traction using consistent displacement jump and nodal reaction forces
*/

#pragma once

#include "PoroCZMComputeLocalTractionBase.h"

class PoroSlipWeakening2dSemiPermeable : public PoroCZMComputeLocalTractionBase
{
public:
  static InputParameters validParams();
  PoroSlipWeakening2dSemiPermeable(const InputParameters & parameters);

protected:
  /// method computing the total traction and its derivatives
  void computeInterfaceTractionAndDerivatives() override;

  Real _T2_o;
  Real _mu_d;
  Real _Dc;

  const MaterialProperty<Real> &_biot_coefficient;
  const MaterialProperty<Real> &_Transmissibility;

  const MaterialProperty<RankTwoTensor> & _rot;
  const VariableValue & _nodal_area;
  const MaterialProperty<Real> & _rhof;
  const MaterialProperty<Real> & _density;

  const MaterialProperty<RankTwoTensor> & _stress;
  const MaterialProperty<RankTwoTensor> & _dstress;

  MaterialProperty<RealVectorValue> & _cumulative_slip;
  const MaterialProperty<RealVectorValue> & _cumulative_slip_old;


  const VariableValue & _interface_pressure_plus; 
  const VariableValue & _interface_pressure_minus;

  const MaterialProperty<RealVectorValue> & _interface_displacement_jump_old;
  const MaterialProperty<RealVectorValue> & _interface_displacement_jump_older;
  
  const VariableValue & _reaction_x;
  const VariableValue & _reaction_neighbor_x;
  const VariableValue & _reaction_y;
  const VariableValue & _reaction_neighbor_y;

  MaterialProperty<Real> & _across_flux_main;
  MaterialProperty<Real> & _across_flux_sec;

};