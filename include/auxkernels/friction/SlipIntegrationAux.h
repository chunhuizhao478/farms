//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "AuxKernel.h"

/**
 * SlipIntegrationAux integrates slip rate to update accumulated slip.
 *
 * S_new = S_old + dt * V
 *
 * This is designed for the Friction SubApp in the MultiApp SEAS architecture.
 * It operates on elements (not boundaries) since the SubApp uses a lower-dimensional
 * fault mesh.
 *
 * Uses coupledValueOld for proper state preservation in MultiApp context.
 * The slip_old_var should be coupled to the slip variable itself.
 */
class SlipIntegrationAux : public AuxKernel
{
public:
  static InputParameters validParams();

  SlipIntegrationAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  /// Slip rate V
  const VariableValue & _slip_rate;

  /// Old value of slip (from previous timestep via coupledValueOld)
  const VariableValue & _slip_old;
};
