#pragma once

#include "InterfaceKernel.h"

class InitialShearTractionKernel : public InterfaceKernel
{
public:
  static InputParameters validParams();

  InitialShearTractionKernel(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual(Moose::DGResidualType type) override;
  const VariableValue & _ini_shear_sts_perturb;
};