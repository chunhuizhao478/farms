#pragma once

#include "InterfaceKernel.h"

/**
 * Rate-and-State Friction Interface Displacement Transfer
 */
class FarmsInterfaceKernelGlobaly : public InterfaceKernel
{
public:
  static InputParameters validParams();

  FarmsInterfaceKernelGlobaly(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual(Moose::DGResidualType type) override;
  //virtual Real computeQpJacobian(Moose::DGJacobianType type) override;
  
  /* the displacements on the plus side (global) */
  const MaterialProperty<RealVectorValue> & _displacements_plus_global;

  /* the displacements on the minus side (global) */
  const MaterialProperty<RealVectorValue> & _displacements_minus_global; 

};