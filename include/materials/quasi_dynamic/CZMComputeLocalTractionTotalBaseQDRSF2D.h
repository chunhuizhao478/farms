#pragma once

#include "CZMComputeLocalTractionBaseQDRSF2D.h"

class CZMComputeLocalTractionTotalBaseQDRSF2D : public CZMComputeLocalTractionBaseQDRSF2D
{
public:
  static InputParameters validParams();
  CZMComputeLocalTractionTotalBaseQDRSF2D(const InputParameters & parameters);
};