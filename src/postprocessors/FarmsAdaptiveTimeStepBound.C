//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FarmsAdaptiveTimeStepBound.h"

registerMooseObject("farmsApp", FarmsAdaptiveTimeStepBound);

InputParameters
FarmsAdaptiveTimeStepBound::validParams()
{
  InputParameters params = GeneralPostprocessor::validParams();

  params.addClassDescription(
      "Returns an adaptive time step bound based on a criterion postprocessor. "
      "Switches between two dt bounds depending on whether the criterion exceeds a threshold. "
      "Useful for adaptive dynamic/quasi-dynamic switching where different modes require "
      "different time step constraints.");

  params.addRequiredParam<PostprocessorName>(
      "criterion_postprocessor",
      "The postprocessor to monitor (e.g., max_dev_strain_rate).");

  params.addRequiredParam<Real>(
      "threshold",
      "The threshold value for switching between time step bounds.");

  params.addRequiredParam<Real>(
      "dt_bound_below_threshold",
      "Time step bound to use when criterion is below/equal threshold (quasi-dynamic mode).");

  params.addRequiredParam<Real>(
      "dt_bound_above_threshold",
      "Time step bound to use when criterion is above threshold (dynamic mode).");

  MooseEnum comparison_type("greater_than less_than greater_equal less_equal", "greater_than");
  params.addParam<MooseEnum>(
      "comparison_type",
      comparison_type,
      "Type of comparison: if criterion > threshold, use dt_bound_above_threshold");

  return params;
}

FarmsAdaptiveTimeStepBound::FarmsAdaptiveTimeStepBound(const InputParameters & parameters)
  : GeneralPostprocessor(parameters),
    _criterion_pp(getPostprocessorValue("criterion_postprocessor")),
    _threshold(getParam<Real>("threshold")),
    _dt_bound_below(getParam<Real>("dt_bound_below_threshold")),
    _dt_bound_above(getParam<Real>("dt_bound_above_threshold")),
    _comparison_type([&]() {
      const auto & comp_str = getParam<MooseEnum>("comparison_type");
      if (comp_str == "greater_than")
        return ComparisonType::GREATER_THAN;
      else if (comp_str == "less_than")
        return ComparisonType::LESS_THAN;
      else if (comp_str == "greater_equal")
        return ComparisonType::GREATER_EQUAL;
      else if (comp_str == "less_equal")
        return ComparisonType::LESS_EQUAL;
      else
        mooseError("Invalid comparison_type");
    }())
{
}

PostprocessorValue
FarmsAdaptiveTimeStepBound::getValue() const
{
  // Determine if condition is met
  bool condition_met = false;

  switch (_comparison_type)
  {
    case ComparisonType::GREATER_THAN:
      condition_met = (_criterion_pp > _threshold);
      break;
    case ComparisonType::LESS_THAN:
      condition_met = (_criterion_pp < _threshold);
      break;
    case ComparisonType::GREATER_EQUAL:
      condition_met = (_criterion_pp >= _threshold);
      break;
    case ComparisonType::LESS_EQUAL:
      condition_met = (_criterion_pp <= _threshold);
      break;
  }

  // Return appropriate bound
  if (condition_met)
    return _dt_bound_above;
  else
    return _dt_bound_below;
}
