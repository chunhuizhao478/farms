//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FirstStepElementIntegralVariablePostprocessor.h"
#include "FEProblem.h"
#include <cmath>

registerMooseObject("farmsApp", FirstStepElementIntegralVariablePostprocessor);

InputParameters
FirstStepElementIntegralVariablePostprocessor::validParams()
{
  InputParameters params = ElementIntegralPostprocessor::validParams();
  params.addRequiredCoupledVar("variable", "The name of the variable that this object operates on");
  params.addClassDescription("Computes a volume integral of the specified variable");
  params.addParam<bool>(
      "use_absolute_value", false, "Whether to use absolute value of the variable or not");
  return params;
}

FirstStepElementIntegralVariablePostprocessor::FirstStepElementIntegralVariablePostprocessor(
    const InputParameters & parameters)
  : ElementIntegralPostprocessor(parameters),
    MooseVariableInterface<Real>(this,
                                 false,
                                 "variable",
                                 Moose::VarKindType::VAR_ANY,
                                 Moose::VarFieldType::VAR_FIELD_STANDARD),
    _u(coupledValue("variable")),
    _u_old(coupledValueOld("variable")),
    _grad_u(coupledGradient("variable")),
    _use_abs_value(getParam<bool>("use_absolute_value")),
    _first_value(0.0),
    _captured(false)
{
  addMooseVariableDependency(&mooseVariableField());
}

Real
FirstStepElementIntegralVariablePostprocessor::computeQpIntegral()
{
  // Always integrate the current field; we'll freeze the first-step value in getValue()
  if (_use_abs_value)
    return std::abs(_u[_qp]);
  else
    return _u[_qp];
}

void
FirstStepElementIntegralVariablePostprocessor::initialize()
{
  // Only perform base initialization before we've captured the first-step integral
  if (!_captured)
    ElementIntegralPostprocessor::initialize();
}

void
FirstStepElementIntegralVariablePostprocessor::execute()
{
  // Only accumulate during the first step
  if (!_captured)
    ElementIntegralPostprocessor::execute();
}

void
FirstStepElementIntegralVariablePostprocessor::finalize()
{
  if (!_captured)
  {
    // Finalize once and cache the first-step integral value
    ElementIntegralPostprocessor::finalize();
    _first_value = ElementIntegralPostprocessor::getValue();
    _captured = true;
  }
  // After capture, do nothing here so the cached value is reused
}

Real
FirstStepElementIntegralVariablePostprocessor::getValue() const
{
  // After capture, always return the cached first-step value
  if (_captured)
    return _first_value;

  // Before capture, return the current base-class value (during first step)
  return ElementIntegralPostprocessor::getValue();
}