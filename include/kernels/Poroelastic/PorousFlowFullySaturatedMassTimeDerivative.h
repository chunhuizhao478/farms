
#pragma once

#include "TimeKernel.h"

class PorousFlowFullySaturatedMassTimeDerivative : public TimeKernel
{
public:
  static InputParameters validParams();

  PorousFlowFullySaturatedMassTimeDerivative(const InputParameters & parameters);

  virtual void computeJacobian() override;

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  bool _lumping;

   /// Constant Biot modulus
  const MaterialProperty<Real> & _biot_modulus;

};