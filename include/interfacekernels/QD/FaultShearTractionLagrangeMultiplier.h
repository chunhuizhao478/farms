#pragma once

#include "InterfaceKernel.h"

/**
 * Interface kernel that enforces shear traction continuity with Lagrange multiplier
 */
class FaultShearTractionLagrangeMultiplier : public InterfaceKernel
{
public:
  static InputParameters validParams();
  FaultShearTractionLagrangeMultiplier(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual(Moose::DGResidualType type) override;
  virtual Real computeQpJacobian(Moose::DGJacobianType type) override;
  virtual Real computeQpOffDiagJacobian(Moose::DGJacobianType type, unsigned int jvar) override;

  /// Lagrange multiplier for shear traction (element side)
  const VariableValue & _lambda_x;
  
  /// Lagrange multiplier for shear traction (neighbor side)
  const VariableValue & _lambda_x_neighbor;
  
  /// Variable number for lambda_x (element side)
  unsigned int _lambda_x_var;
  
  /// Variable number for lambda_x (neighbor side)
  unsigned int _lambda_x_neighbor_var;
};