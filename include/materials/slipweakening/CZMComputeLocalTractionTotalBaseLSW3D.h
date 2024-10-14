#pragma once

#include "CZMComputeLocalTractionBaseLSW3D.h"

class CZMComputeLocalTractionTotalBaseLSW3D : public CZMComputeLocalTractionBaseLSW3D
{
public:
  static InputParameters validParams();
  CZMComputeLocalTractionTotalBaseLSW3D(const InputParameters & parameters);
};