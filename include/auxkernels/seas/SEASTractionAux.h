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
 * SEASTractionAux reads the elastic traction from an InterfaceMaterial
 * and stores it in an AuxVariable for use by SEASSlipRateAux.
 *
 * This works on boundary sidesets where an InterfaceMaterial provides
 * the traction property.
 *
 * This is part of the staggered SEAS solver approach.
 */
class SEASTractionAux : public AuxKernel
{
public:
  static InputParameters validParams();

  SEASTractionAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  /// Traction from interface material
  const MaterialProperty<Real> & _traction_prop;

  /// Pre-stress (background traction)
  const Real _tau_pre;
};
