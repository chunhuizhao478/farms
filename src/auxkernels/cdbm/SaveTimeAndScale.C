#include "SaveTimeAndScale.h"

registerMooseObject("farmsApp", SaveTimeAndScale);

InputParameters
SaveTimeAndScale::validParams()
{
  InputParameters params = AuxKernel::validParams();
  return params;
}

SaveTimeAndScale::SaveTimeAndScale(const InputParameters & parameters)
  : AuxKernel(parameters)

{
}

Real
SaveTimeAndScale::computeValue()
{

  Real scale = 1e11;
  return _t / scale;
}