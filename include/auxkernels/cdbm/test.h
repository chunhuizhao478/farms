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
 * Constant auxiliary value
 */
class test : public AuxKernel
{
public:
  static InputParameters validParams();

  test(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  /// The value being set for the current node/element
  const Real & _value;
};