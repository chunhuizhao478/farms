//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#pragma once

#include "DegradationFunctionBase.h"

/**
 * Predefined rational degradation function
 */
class RationalDegradationFunctionCZM : public DegradationFunctionBase
{
public:
  static InputParameters validParams();

  RationalDegradationFunctionCZM(const InputParameters & parameters);
};
