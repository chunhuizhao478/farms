#pragma once

#include "CZMComputeLocalTractionBaseRSF2D.h"

class CZMComputeLocalTractionTotalBaseRSF2D : public CZMComputeLocalTractionBaseRSF2D
{
public:
  static InputParameters validParams();
  CZMComputeLocalTractionTotalBaseRSF2D(const InputParameters & parameters);
};