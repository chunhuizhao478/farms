#pragma once

#include "InterfaceKernel.h"

/**
 * Rate-and-State Friction Interface Displacement Transfer
 */
class RateStateInterfaceKernelGlobalx : public InterfaceKernel
{
public:
  static InputParameters validParams();

  RateStateInterfaceKernelGlobalx(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual(Moose::DGResidualType type) override;
  //virtual Real computeQpJacobian(Moose::DGJacobianType type) override;
  
  const MaterialProperty<Real> & _disp_strike_plus;
  const MaterialProperty<Real> & _disp_strike_minus;

};