#include "Assembly.h"
#include "CZMComputeLocalTractionTotalBaseQDRSF2D.h"

InputParameters
CZMComputeLocalTractionTotalBaseQDRSF2D::validParams()
{
  InputParameters params = CZMComputeLocalTractionBaseQDRSF2D::validParams();

  params.addClassDescription("Base class for implementing cohesive zone constituive material "
                             "models that can be formulated using the total displacement jump");
  return params;
}

CZMComputeLocalTractionTotalBaseQDRSF2D::CZMComputeLocalTractionTotalBaseQDRSF2D(
    const InputParameters & parameters)
  : CZMComputeLocalTractionBaseQDRSF2D(parameters)
{
}