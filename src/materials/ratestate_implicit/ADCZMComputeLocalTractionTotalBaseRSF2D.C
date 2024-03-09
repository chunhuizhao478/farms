//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADCZMComputeLocalTractionTotalBaseRSF2D.h"
//#include "CZMComputeLocalTractionTotalBase.h"

InputParameters
ADCZMComputeLocalTractionTotalBaseRSF2D::validParams()
{
  InputParameters params = ADCZMComputeLocalTractionTotalBaseRSF2D::validParams();
  return params;
}

ADCZMComputeLocalTractionTotalBaseRSF2D::ADCZMComputeLocalTractionTotalBaseRSF2D(
    const InputParameters & parameters)
  : ADCZMComputeLocalTractionBaseRSF2D(parameters)
{
}