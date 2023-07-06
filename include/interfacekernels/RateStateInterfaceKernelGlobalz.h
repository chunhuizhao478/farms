#pragma once

#include "InterfaceKernel.h"

/**
 * Rate-and-State Friction Interface Displacement Transfer
 */
class RateStateInterfaceKernelGlobalz : public InterfaceKernel
{
public:
  static InputParameters validParams();

  RateStateInterfaceKernelGlobalz(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual(Moose::DGResidualType type) override;
  //virtual Real computeQpJacobian(Moose::DGJacobianType type) override;
  
  const MaterialProperty<Real> & _disp_dip_plus;
  const MaterialProperty<Real> & _disp_dip_minus;

};