#pragma once

#include "PoroCZMComputeLocalTractionBaseRSF2D.h"

class PoroCZMComputeLocalTractionTotalBaseRSF2D : public PoroCZMComputeLocalTractionBaseRSF2D
{
public:
  static InputParameters validParams();
  PoroCZMComputeLocalTractionTotalBaseRSF2D(const InputParameters & parameters);
};