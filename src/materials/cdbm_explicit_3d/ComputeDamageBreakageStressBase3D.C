//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ComputeDamageBreakageStressBase3D.h"
#include "ComputeElasticityTensorBase.h"
#include "Function.h"

InputParameters
ComputeDamageBreakageStressBase3D::validParams()
{
  InputParameters params = ComputeGeneralDamageBreakageStressBase3D::validParams();
  params.suppressParameter<bool>("use_displaced_mesh");
  return params;
}

ComputeDamageBreakageStressBase3D::ComputeDamageBreakageStressBase3D(const InputParameters & parameters)
  : ComputeGeneralDamageBreakageStressBase3D(parameters)
{
  if (getParam<bool>("use_displaced_mesh"))
    mooseError("The stress calculator needs to run on the undisplaced mesh.");
}
