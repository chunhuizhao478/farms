#include "Assembly.h"
#include "IMComputeLocalTractionTotalBaseRSF2D.h"

InputParameters
IMComputeLocalTractionTotalBaseRSF2D::validParams()
{
  InputParameters params = IMComputeLocalTractionBaseRSF2D::validParams();

  params.addClassDescription("Base class for implementing cohesive zone constituive material "
                             "models that can be formulated using the total displacement jump");
  return params;
}

IMComputeLocalTractionTotalBaseRSF2D::IMComputeLocalTractionTotalBaseRSF2D(
    const InputParameters & parameters)
  : IMComputeLocalTractionBaseRSF2D(parameters)
{
}