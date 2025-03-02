#pragma once

#include "CZMComputeLocalTractionBaseSWF2D.h"

class CZMComputeLocalTractionTotalBaseSWF2D : public CZMComputeLocalTractionBaseSWF2D
{
public:
  static InputParameters validParams();
  CZMComputeLocalTractionTotalBaseSWF2D(const InputParameters & parameters);
};