#include "Assembly.h"
#include "CZMComputeLocalTractionTotalBaseRSF3D.h"

InputParameters
CZMComputeLocalTractionTotalBaseRSF3D::validParams()
{
  InputParameters params = CZMComputeLocalTractionBaseRSF3D::validParams();

  params.addClassDescription("Base class for implementing cohesive zone constituive material "
                             "models that can be formulated using the total displacement jump");
  return params;
}

CZMComputeLocalTractionTotalBaseRSF3D::CZMComputeLocalTractionTotalBaseRSF3D(
    const InputParameters & parameters)
  : CZMComputeLocalTractionBaseRSF3D(parameters)
{
}