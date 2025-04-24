//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "test.h"

registerMooseObject("farmsApp", test);

InputParameters
test::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("Creates a constant field in the domain.");
  params.addParam<Real>("value", 0.0, "Some constant value that can be read from the input file");
  params.declareControllable("value");
  return params;
}

test::test(const InputParameters & parameters)
  : AuxKernel(parameters), _value(getParam<Real>("value"))
{
}

Real
test::computeValue()
{

  // Get coordinates (no rotation applied)
  Real xcoord = _q_point[_qp](0); // strike direction
  Real ycoord = _q_point[_qp](1); // normal direction

  std::cout << "xcoord: " << xcoord << std::endl;
  std::cout << "ycoord: " << ycoord << std::endl;

  return _value;
}