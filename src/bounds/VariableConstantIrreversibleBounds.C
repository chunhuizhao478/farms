//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "VariableConstantIrreversibleBounds.h"

registerMooseObject("farmsApp", VariableConstantIrreversibleBounds);

InputParameters
VariableConstantIrreversibleBounds::validParams()
{
  InputParameters params = VariableBoundsBase::validParams();
  params.addClassDescription(
      "Provides variable bound of a variable for the PETSc's variational inequalities solver initially, and a old value bound after the first step to ensure irreversibility.");
  params.addRequiredCoupledVar("bound_value", "The value of bound for the variable");
  return params;
}

VariableConstantIrreversibleBounds::VariableConstantIrreversibleBounds(const InputParameters & parameters)
  : VariableBoundsBase(parameters), _bound_value(coupledValue("bound_value"))
{
}