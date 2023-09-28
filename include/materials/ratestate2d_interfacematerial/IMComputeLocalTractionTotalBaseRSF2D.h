#pragma once

#include "IMComputeLocalTractionBaseRSF2D.h"

class IMComputeLocalTractionTotalBaseRSF2D : public IMComputeLocalTractionBaseRSF2D
{
public:
  static InputParameters validParams();
  IMComputeLocalTractionTotalBaseRSF2D(const InputParameters & parameters);
};