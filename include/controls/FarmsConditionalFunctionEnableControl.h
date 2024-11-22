//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "FarmsConditionalEnableControl.h"

class Function;

/**
 * Control for enabling/disabling objects when a function value is true
 */
class FarmsConditionalFunctionEnableControl : public FarmsConditionalEnableControl
{
public:
  static InputParameters validParams();

  FarmsConditionalFunctionEnableControl(const InputParameters & parameters);

protected:
  virtual bool conditionMet(const unsigned int & i) override;

private:
  /// The function to give a true or false value
  const Function & _function;
};