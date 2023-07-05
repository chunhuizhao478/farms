#include "Assembly.h"
#include "CZMComputeLocalTractionTotalBaseRSF2D.h"

InputParameters
CZMComputeLocalTractionTotalBaseRSF2D::validParams()
{
  InputParameters params = CZMComputeLocalTractionBaseRSF2D::validParams();

  params.addClassDescription("Base class for implementing cohesive zone constituive material "
                             "models that can be formulated using the total displacement jump");
  return params;
}

CZMComputeLocalTractionTotalBaseRSF2D::CZMComputeLocalTractionTotalBaseRSF2D(
    const InputParameters & parameters)
  : CZMComputeLocalTractionBaseRSF2D(parameters)
{
}