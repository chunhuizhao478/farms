#pragma once

#include "CZMComputeLocalTractionTotalBase.h"

/**
 * Implementation of the Slip Weakening Law
 **/
class SlipWeakeningFriction2dxx : public CZMComputeLocalTractionTotalBase
{
public:
  static InputParameters validParams();
  SlipWeakeningFriction2dxx(const InputParameters & parameters);

protected:
  /// method computing the total traction and its derivatives
  void computeInterfaceTractionAndDerivatives() override;

  Real _T2_o;
  const VariableValue & _nodal_area;
  Real _mu_d;
  Real _Dc;

  const MaterialProperty<RankTwoTensor> & _stress;
  const MaterialProperty<RankTwoTensor> & _dstress;
  const MaterialProperty<RealVectorValue> & _interface_displacement_jump_old;
  const MaterialProperty<RealVectorValue> & _interface_displacement_jump_older;
  const MaterialProperty<Real> & _density;


};