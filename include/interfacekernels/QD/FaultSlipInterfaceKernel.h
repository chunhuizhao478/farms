//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "InterfaceKernel.h"

/**
 * Interface kernel that enforces a slip condition between
 * two surfaces of a fault with penalty enforcement
 */
class FaultSlipInterfaceKernel : public InterfaceKernel
{
public:
  static InputParameters validParams();

  FaultSlipInterfaceKernel(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual(Moose::DGResidualType type) override;
  virtual Real computeQpJacobian(Moose::DGJacobianType type) override;
  virtual Real computeQpOffDiagJacobian(Moose::DGJacobianType type, unsigned int jvar) override;

  /// Coupled slip value that defines the target jump across the interface
  const VariableValue & _slip;

  /// Variable number for the slip
  const unsigned int _slip_var;

  /// Coupled displacement on the element side
  const VariableValue & _coupled_disp;

  /// Variable number for the coupled displacement
  const unsigned int _coupled_disp_num;

  /// Coupled displacement on the neighbor side
  const VariableValue & _coupled_disp_neighbor;

  /// Variable number for the coupled displacement on the neighbor side
  const unsigned int _coupled_disp_neighbor_num;
};