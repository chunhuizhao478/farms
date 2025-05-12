#include "NDErrorPsiMeasure.h"

registerMooseObject("farmsApp", NDErrorPsiMeasure);

InputParameters
NDErrorPsiMeasure::validParams()
{
  InputParameters params = AuxKernel::validParams();
  return params;
}

NDErrorPsiMeasure::NDErrorPsiMeasure(const InputParameters & parameters)
  : AuxKernel(parameters),
  _psie_active(getMaterialProperty<Real>("psie_active")),
  _psie_active_old(getMaterialPropertyOldByName<Real>("psie_active"))
{
}

Real
NDErrorPsiMeasure::computeValue()
{
  return _psie_active[_qp] - _psie_active_old[_qp];
}