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
class PoroSlipWeakening2dImpermeable2 : public CZMComputeLocalTractionTotalBase
{
public:
  static InputParameters validParams();
  PoroSlipWeakening2dImpermeable2(const InputParameters & parameters);

protected:
  /// method computing the total traction and its derivatives
  void computeInterfaceTractionAndDerivatives() override;

  Real _T2_o;
  Real _mu_d;
  Real _Dc;

  const VariableValue & _interface_pressure_plus; 
  const VariableValue & _interface_pressure_minus;

  const MaterialProperty<RankTwoTensor> & _stress;
  const MaterialProperty<RankTwoTensor> & _dstress;
  const MaterialProperty<RealVectorValue> & _interface_displacement_jump_old;
  const MaterialProperty<RealVectorValue> & _interface_displacement_jump_older;

  const VariableValue & _nodal_area;
  const MaterialProperty<Real> & _density;
  

};
