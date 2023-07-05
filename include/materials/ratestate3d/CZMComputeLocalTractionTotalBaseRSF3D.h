#pragma once

#include "CZMComputeLocalTractionBaseRSF3D.h"

class CZMComputeLocalTractionTotalBaseRSF3D : public CZMComputeLocalTractionBaseRSF3D
{
public:
  static InputParameters validParams();
  CZMComputeLocalTractionTotalBaseRSF3D(const InputParameters & parameters);
};