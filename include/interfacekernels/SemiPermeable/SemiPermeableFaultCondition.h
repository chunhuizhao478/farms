#pragma once

#include "InterfaceKernel.h"


class SemiPermeableFaultCondition : public InterfaceKernel
{
public:
  static InputParameters validParams();

  SemiPermeableFaultCondition(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual(Moose::DGResidualType type) override;
  virtual Real computeQpJacobian(Moose::DGJacobianType type) override;

  const MaterialProperty<Real> & _trans;

  const VariableValue & _pore_pressure_main; 
  const VariableValue & _pore_pressure_secondary;

};