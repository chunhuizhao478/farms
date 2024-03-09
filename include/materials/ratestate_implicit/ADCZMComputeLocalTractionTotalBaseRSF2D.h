//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ADCZMComputeLocalTractionBaseRSF2D.h"

/**
 * AD equivalent of CZMComputeLocalTractionTotalBase
 */
class ADCZMComputeLocalTractionTotalBaseRSF2D : public ADCZMComputeLocalTractionBaseRSF2D
{
public:
  static InputParameters validParams();
  ADCZMComputeLocalTractionTotalBaseRSF2D(const InputParameters & parameters);
};