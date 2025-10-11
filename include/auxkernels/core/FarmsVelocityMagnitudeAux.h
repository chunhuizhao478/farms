#pragma once

#include "AuxKernel.h"

/**
 * Computes the magnitude of the velocity vector from its components.
 */
class FarmsVelocityMagnitudeAux : public AuxKernel
{
public:
  static InputParameters validParams();

  FarmsVelocityMagnitudeAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  const VariableValue & _vel_x;
  const VariableValue & _vel_y;
  const VariableValue * _vel_z;
  const bool _has_z;
};
