#include "ErrorPsiMeasure.h"

registerMooseObject("farmsApp", ErrorPsiMeasure);

InputParameters
ErrorPsiMeasure::validParams()
{
  InputParameters params = AuxKernel::validParams();
  return params;
}

ErrorPsiMeasure::ErrorPsiMeasure(const InputParameters & parameters)
  : AuxKernel(parameters),
  _psie_active(getADMaterialProperty<Real>("psie_active")),
  _psie_active_old(getMaterialPropertyOldByName<Real>("psie_active"))
{
}

Real
ErrorPsiMeasure::computeValue()
{
  return MetaPhysicL::raw_value(_psie_active[_qp]) - _psie_active_old[_qp];
}