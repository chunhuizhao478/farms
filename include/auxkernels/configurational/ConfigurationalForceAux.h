//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "AuxKernel.h"

// Forward declaration
class ConfigurationalForceUserObject;

/**
 * AuxKernel to output configurational force components computed by
 * ConfigurationalForceUserObject to AuxVariables for visualization.
 */
class ConfigurationalForceAux : public AuxKernel
{
public:
  static InputParameters validParams();

  ConfigurationalForceAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  /// Reference to the UserObject that computed the configurational forces
  const ConfigurationalForceUserObject & _config_force_uo;

  /// Which component to output (0=x, 1=y, 2=z, 3=magnitude)
  const unsigned int _component;
};
