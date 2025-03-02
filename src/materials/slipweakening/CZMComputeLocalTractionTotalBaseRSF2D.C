#include "Assembly.h"
#include "CZMComputeLocalTractionTotalBaseSWF2D.h"

InputParameters
CZMComputeLocalTractionTotalBaseSWF2D::validParams()
{
  InputParameters params = CZMComputeLocalTractionBaseSWF2D::validParams();

  params.addClassDescription("Base class for implementing cohesive zone constituive material "
                             "models that can be formulated using the total displacement jump");
  return params;
}

CZMComputeLocalTractionTotalBaseSWF2D::CZMComputeLocalTractionTotalBaseSWF2D(
    const InputParameters & parameters)
  : CZMComputeLocalTractionBaseSWF2D(parameters)
{
}