#include "Assembly.h"
#include "PoroCZMComputeLocalTractionTotalBaseRSF2D.h"

InputParameters
PoroCZMComputeLocalTractionTotalBaseRSF2D::validParams()
{
  InputParameters params = PoroCZMComputeLocalTractionBaseRSF2D::validParams();

  params.addClassDescription("Base class for implementing cohesive zone constituive material "
                             "models that can be formulated using the total displacement jump");
  return params;
}

PoroCZMComputeLocalTractionTotalBaseRSF2D::PoroCZMComputeLocalTractionTotalBaseRSF2D(
    const InputParameters & parameters)
  : PoroCZMComputeLocalTractionBaseRSF2D(parameters)
{
}