#pragma once

#include "InterfaceKernel.h"

/**
 * Rate-and-State Friction Interface Displacement Transfer
 */
class RateStateInterfaceKernelGlobalxdev : public InterfaceKernel
{
public:
  static InputParameters validParams();

  RateStateInterfaceKernelGlobalxdev(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual(Moose::DGResidualType type) override;
  //virtual Real computeQpJacobian(Moose::DGJacobianType type) override;

  const VariableValue & _disp_strike_plus;
  const VariableValue & _disp_strike_minus;

};