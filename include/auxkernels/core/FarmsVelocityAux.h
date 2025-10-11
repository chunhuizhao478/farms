#pragma once

#include "AuxKernel.h"

/**
 * Computes the time derivative of a displacement variable and stores it in an auxiliary velocity
 * variable. Intended for use in quasi-dynamic simulations where the velocity field is needed but
 * inertia is neglected.
 */
class FarmsVelocityAux : public AuxKernel
{
public:
  static InputParameters validParams();

  FarmsVelocityAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  /// Time derivative of the coupled displacement variable
  const VariableValue & _disp_dot;
};
