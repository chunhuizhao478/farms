#pragma once

#include "CZMComputeLocalTractionTotalBase.h"

/**
 * Implementation of the Slip Weakening Law
 **/
class SlipWeakeningFriction3d : public CZMComputeLocalTractionTotalBase
{
public:
  static InputParameters validParams();
  SlipWeakeningFriction3d(const InputParameters & parameters);

protected:
  /// method computing the total traction and its derivatives
  void computeInterfaceTractionAndDerivatives() override;

  Real _T2_o;
  Real _area;
  Real _mu_d;
  Real _Dc;

  const MaterialProperty<RankTwoTensor> & _stress;
  const MaterialProperty<RealVectorValue> & _interface_displacement_jump_old;
  const MaterialProperty<RealVectorValue> & _interface_displacement_jump_older;
  const MaterialProperty<Real> & _density;


};