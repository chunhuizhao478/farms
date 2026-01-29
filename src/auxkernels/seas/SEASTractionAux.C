//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SEASTractionAux.h"

registerMooseObject("farmsApp", SEASTractionAux);

InputParameters
SEASTractionAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription(
      "Reads elastic traction from an InterfaceMaterial property. "
      "Part of staggered SEAS solver.");

  params.addRequiredParam<MaterialPropertyName>("traction_property",
                                                 "Name of the traction material property");
  params.addParam<Real>("tau_pre", 0.0, "Pre-stress / background traction (Pa)");

  return params;
}

SEASTractionAux::SEASTractionAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _traction_prop(getMaterialProperty<Real>("traction_property")),
    _tau_pre(getParam<Real>("tau_pre"))
{
}

Real
SEASTractionAux::computeValue()
{
  return _traction_prop[_qp] + _tau_pre;
}
