
#pragma once

#include "TimeKernel.h"

class PorousFlowFullySaturatedVolumetricExpansion : public TimeKernel
{
public:
  static InputParameters validParams();

  PorousFlowFullySaturatedVolumetricExpansion(const InputParameters & parameters);

  virtual void computeJacobian() override;

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  bool _lumping;

  /// Biot coefficient (used in simulations involving Mechanical deformations)
  const Real _biot_coefficient;
};