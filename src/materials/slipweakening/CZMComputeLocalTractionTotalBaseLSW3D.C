#include "Assembly.h"
#include "CZMComputeLocalTractionTotalBaseLSW3D.h"

InputParameters
CZMComputeLocalTractionTotalBaseLSW3D::validParams()
{
  InputParameters params = CZMComputeLocalTractionBaseLSW3D::validParams();

  params.addClassDescription("Base class for implementing cohesive zone constituive material "
                             "models that can be formulated using the total displacement jump");
  return params;
}

CZMComputeLocalTractionTotalBaseLSW3D::CZMComputeLocalTractionTotalBaseLSW3D(
    const InputParameters & parameters)
  : CZMComputeLocalTractionBaseLSW3D(parameters)
{
}