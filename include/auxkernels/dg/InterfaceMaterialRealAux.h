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
 * InterfaceMaterialRealAux reads a Real material property from an interface
 * and stores it in an auxiliary variable on the boundary.
 *
 * This is useful for outputting interface material properties (e.g., slip,
 * slip rate, state variable) computed by InterfaceMaterial objects like
 * DGRateStateFrictionMaterial.
 */
class InterfaceMaterialRealAux : public AuxKernel
{
public:
  static InputParameters validParams();

  InterfaceMaterialRealAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  /// The interface material property to read
  const MaterialProperty<Real> & _prop;
};
