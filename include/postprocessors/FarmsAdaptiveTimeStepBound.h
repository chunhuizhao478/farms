//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "GeneralPostprocessor.h"

/**
 * Returns an adaptive time step bound based on a criterion postprocessor value.
 *
 * This postprocessor switches between two different time step bounds based on
 * whether a monitored postprocessor (e.g., deviatoric strain rate) exceeds a threshold.
 *
 * Use case: Return large dt_max for quasi-dynamic mode, small dt_max for dynamic mode.
 */
class FarmsAdaptiveTimeStepBound : public GeneralPostprocessor
{
public:
  static InputParameters validParams();

  FarmsAdaptiveTimeStepBound(const InputParameters & parameters);

  virtual void initialize() override {}
  virtual void execute() override {}
  virtual PostprocessorValue getValue() const override;

protected:
  /// The postprocessor to monitor (e.g., max deviatoric strain rate)
  const PostprocessorValue & _criterion_pp;

  /// Threshold value for switching
  const Real _threshold;

  /// Time step bound when criterion is below threshold (quasi-dynamic mode)
  const Real _dt_bound_below;

  /// Time step bound when criterion is above threshold (dynamic mode)
  const Real _dt_bound_above;

  /// Comparison type
  enum class ComparisonType
  {
    GREATER_THAN,
    LESS_THAN,
    GREATER_EQUAL,
    LESS_EQUAL
  };

  const ComparisonType _comparison_type;
};
