//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MeshSize.h"
#include <cmath>
#include <random>
#include <algorithm>

registerMooseObject("farmsApp", MeshSize);

InputParameters
MeshSize::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("AuxKernel that computes element size");
  return params;
}

MeshSize::MeshSize(const InputParameters & parameters) : 
AuxKernel(parameters)
{
}

Real
MeshSize::computeValue()
{
  //compute element length assume equal length tetrahedron
  Real length_a = pow(6*sqrt(2)*_current_elem_volume,1.0/3.0);

  return length_a;
}