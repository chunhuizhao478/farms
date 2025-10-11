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

/**
 * Control that enables/disables objects based on a postprocessor value comparison
 *
 * This control evaluates a condition comparing a postprocessor value against a threshold,
 * and enables/disables specified objects accordingly. Supports bidirectional switching
 * via the reverse_on_false parameter.
 */
class FarmsConditionalPostprocessorEnableControl : public FarmsConditionalEnableControl
{
public:
  static InputParameters validParams();

  FarmsConditionalPostprocessorEnableControl(const InputParameters & parameters);

protected:
  virtual bool conditionMet(const unsigned int & i) override;

  /// Comparison type enumeration
  enum class ComparisonType
  {
    GREATER_THAN,
    LESS_THAN,
    GREATER_EQUAL,
    LESS_EQUAL,
    EQUAL,
    NOT_EQUAL
  };

  /// The postprocessor to monitor
  const PostprocessorValue & _pp_value;

  /// Threshold value for comparison
  const Real _threshold;

  /// Type of comparison to perform
  const ComparisonType _comparison_type;
};
