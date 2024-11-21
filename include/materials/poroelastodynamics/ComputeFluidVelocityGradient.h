
#pragma once

#include "DarcyFluidVelocity.h"

class ComputeFluidVelocityGradient : public DarcyFluidVelocity
{
public:
  static InputParameters validParams();

  ComputeFluidVelocityGradient(const InputParameters & parameters);

  virtual void computeProperties() override;
};