/* 
Material Description of Slip Weakening Friction v2
1. Release the zero normal traction perturbation restriction by passing reaction forces plus/mins surfaces
2. Release the fix nodal area/length value along the fault by using NodalArea (Contact Module) evaluated at the interface
*/

#pragma once

#include "CZMComputeLocalTractionTotalBase.h"

/**
 * Implementation of the Slip Weakening Law
 **/
class SlipWeakeningFriction2dv2 : public CZMComputeLocalTractionTotalBase
{
public:
  static InputParameters validParams();
  SlipWeakeningFriction2dv2(const InputParameters & parameters);

protected:
  /// method computing the total traction and its derivatives
  void computeInterfaceTractionAndDerivatives() override;

  Real _T2_o;
  Real _mu_d;
  Real _Dc;

  const MaterialProperty<RankTwoTensor> & _stress;
  const MaterialProperty<RealVectorValue> & _interface_displacement_jump_old;
  const MaterialProperty<RealVectorValue> & _interface_displacement_jump_older;

  const MaterialProperty<Real> & _density;

  const MaterialProperty<RankTwoTensor> & _rot;

  const VariableValue & _reaction_x;
  const VariableValue & _reaction_neighbor_x;
  const VariableValue & _reaction_y;
  const VariableValue & _reaction_neighbor_y;
  const VariableValue & _nodal_area;
 

};
