//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FarmsConditionalPostprocessorEnableControl.h"

registerMooseObject("farmsApp", FarmsConditionalPostprocessorEnableControl);

InputParameters
FarmsConditionalPostprocessorEnableControl::validParams()
{
  InputParameters params = FarmsConditionalEnableControl::validParams();

  params.addClassDescription(
      "Control that enables/disables objects based on a postprocessor value comparison. "
      "Supports bidirectional switching when reverse_on_false = true.");

  params.addRequiredParam<PostprocessorName>(
      "postprocessor",
      "The postprocessor whose value will be compared against the threshold.");

  params.addRequiredParam<Real>(
      "threshold",
      "The threshold value for comparison.");

  MooseEnum comparison_type("greater_than less_than greater_equal less_equal equal not_equal",
                            "greater_than");
  params.addParam<MooseEnum>(
      "comparison_type",
      comparison_type,
      "Type of comparison to perform: greater_than (>), less_than (<), "
      "greater_equal (>=), less_equal (<=), equal (==), not_equal (!=)");

  return params;
}

FarmsConditionalPostprocessorEnableControl::FarmsConditionalPostprocessorEnableControl(
    const InputParameters & parameters)
  : FarmsConditionalEnableControl(parameters),
    _pp_value(getPostprocessorValue("postprocessor")),
    _threshold(getParam<Real>("threshold")),
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
      else if (comp_str == "equal")
        return ComparisonType::EQUAL;
      else if (comp_str == "not_equal")
        return ComparisonType::NOT_EQUAL;
      else
        mooseError("Invalid comparison_type");
    }())
{
}

bool
FarmsConditionalPostprocessorEnableControl::conditionMet(const unsigned int & /* i */)
{
  // Perform comparison based on type
  switch (_comparison_type)
  {
    case ComparisonType::GREATER_THAN:
      return _pp_value > _threshold;
    case ComparisonType::LESS_THAN:
      return _pp_value < _threshold;
    case ComparisonType::GREATER_EQUAL:
      return _pp_value >= _threshold;
    case ComparisonType::LESS_EQUAL:
      return _pp_value <= _threshold;
    case ComparisonType::EQUAL:
      return std::abs(_pp_value - _threshold) < 1e-10;
    case ComparisonType::NOT_EQUAL:
      return std::abs(_pp_value - _threshold) >= 1e-10;
    default:
      mooseError("Unknown comparison type");
  }
}
