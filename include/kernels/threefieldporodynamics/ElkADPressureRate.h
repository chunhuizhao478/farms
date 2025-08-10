//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

/**
 * Created by Chunhui Zhao, Jul 3rd, 2024
 * Implement dot(p) / M term
 * dot(p): pressure time derivative (pressure is the primary variable)
 * M : biot modulus 
 */
#pragma once

#include "ADTimeKernelValue.h"

class ElkADPressureRate : public ADTimeKernelValue
{
public:
  static InputParameters validParams();

  ElkADPressureRate(const InputParameters & parameters);

protected:
  virtual ADReal precomputeQpResidual() override;
};