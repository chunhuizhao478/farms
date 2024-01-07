/*
AuxKernel of Passing Variable
*/

#include "CompArctanyx.h"
#include <cmath>

registerMooseObject("farmsApp", CompArctanyx);

InputParameters
CompArctanyx::validParams()
{
  InputParameters params = AuxKernel::validParams();


  return params;
}

CompArctanyx::CompArctanyx(const InputParameters & parameters)
  : AuxKernel(parameters)

{
}

Real
CompArctanyx::computeValue()
{
  Real xcoord = (*_current_node)(0);
  Real ycoord = (*_current_node)(1);
  Real res = std::atan2(ycoord,xcoord);
  return res;
}