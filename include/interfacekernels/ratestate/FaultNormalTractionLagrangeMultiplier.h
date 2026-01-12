#pragma once

#include "InterfaceKernel.h"

/**
 * Interface kernel that enforces normal traction continuity with Lagrange multiplier
 */
class FaultNormalTractionLagrangeMultiplier : public InterfaceKernel
{
public:
  static InputParameters validParams();
  FaultNormalTractionLagrangeMultiplier(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual(Moose::DGResidualType type) override;
  virtual Real computeQpJacobian(Moose::DGJacobianType type) override;
  virtual Real computeQpOffDiagJacobian(Moose::DGJacobianType type, unsigned int jvar) override;

  /// Lagrange multiplier for normal traction (element side)
  const VariableValue & _lambda_y;
  
  /// Lagrange multiplier for normal traction (neighbor side)
  const VariableValue & _lambda_y_neighbor;
  
  /// Variable number for lambda_y (element side)
  unsigned int _lambda_y_var;
  
  /// Variable number for lambda_y (neighbor side)
  unsigned int _lambda_y_neighbor_var;
};