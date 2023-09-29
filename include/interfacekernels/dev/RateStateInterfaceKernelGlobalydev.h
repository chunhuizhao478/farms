#pragma once

#include "InterfaceKernel.h"

/**
 * Rate-and-State Friction Interface Displacement Transfer
 */
class RateStateInterfaceKernelGlobalydev : public InterfaceKernel
{
public:
  static InputParameters validParams();

  RateStateInterfaceKernelGlobalydev(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual(Moose::DGResidualType type) override;
  //virtual Real computeQpJacobian(Moose::DGJacobianType type) override;
  
  const VariableValue & _disp_normal_plus;
  const VariableValue & _disp_normal_minus;

};