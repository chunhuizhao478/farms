//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "VariableConstantBounds.h"

registerMooseObject("farmsApp", VariableConstantBounds);

InputParameters
VariableConstantBounds::validParams()
{
  InputParameters params = VariableBoundsBase::validParams();
  params.addClassDescription(
      "Provides constant bound of a variable for the PETSc's variational inequalities solver");
  params.addRequiredCoupledVar("bound_value", "The value of bound for the variable");
  return params;
}

VariableConstantBounds::VariableConstantBounds(const InputParameters & parameters)
  : VariableBoundsBase(parameters), _bound_value(coupledValue("bound_value"))
{
}