//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "InterfaceMaterialRealAux.h"

registerMooseObject("farmsApp", InterfaceMaterialRealAux);

InputParameters
InterfaceMaterialRealAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription(
      "Reads a Real material property from an interface material and stores it "
      "in an auxiliary variable. Useful for outputting interface properties like "
      "slip, slip rate, and state variable from DGRateStateFrictionMaterial.");

  params.addRequiredParam<MaterialPropertyName>("property",
                                                "The name of the interface material property");

  return params;
}

InterfaceMaterialRealAux::InterfaceMaterialRealAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _prop(getMaterialProperty<Real>("property"))
{
  // This AuxKernel should be used on boundaries
  if (!isNodal())
  {
    // Element-based aux variable is fine for MONOMIAL
  }
}

Real
InterfaceMaterialRealAux::computeValue()
{
  return _prop[_qp];
}
