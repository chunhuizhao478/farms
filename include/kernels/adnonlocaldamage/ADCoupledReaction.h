#pragma once

#include "ADKernel.h"

/**
 * This kernel implements a simple reaction term with a coupled local equivalent strain using AD
 */
class ADCoupledReaction : public ADKernel
{
public:
  static InputParameters validParams();
  ADCoupledReaction(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  /// Reaction rate
  const Real & _rate;
  
  /// Coupled local equivalent strain
  const ADMaterialProperty<Real> & _eqstrain_local;
  /// Length scale for gradient activity parameter
  const Real _length_scale;
  /// Equivalent strain at which the gradient activity starts
  const Real _kappa_i;  
  /// Minimum value of the gradient activity parameter for the equivalent strain
  const Real _c0;
};