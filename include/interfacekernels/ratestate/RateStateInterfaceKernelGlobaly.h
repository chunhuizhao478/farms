#pragma once

#include "InterfaceKernel.h"

/**
 * Rate-and-State Friction Interface Displacement Transfer
 */
class RateStateInterfaceKernelGlobaly : public InterfaceKernel
{
public:
  static InputParameters validParams();

  RateStateInterfaceKernelGlobaly(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual(Moose::DGResidualType type) override;
  //virtual Real computeQpJacobian(Moose::DGJacobianType type) override;
  
  const MaterialProperty<Real> & _disp_normal_plus;
  const MaterialProperty<Real> & _disp_normal_minus;

};