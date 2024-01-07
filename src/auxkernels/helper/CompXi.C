/*
AuxKernel of Passing Variable
*/

#include "CompXi.h"

registerMooseObject("farmsApp", CompXi);

InputParameters
CompXi::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addRequiredCoupledVar("strainxx","strain xx");
  params.addRequiredCoupledVar("strainxy","strain xy");
  params.addRequiredCoupledVar("strainyy","strain yy");
  return params;
}

CompXi::CompXi(const InputParameters & parameters)
  : AuxKernel(parameters),
  _strainxx(coupledValue("strainxx")),
  _strainxy(coupledValue("strainxy")),
  _strainyy(coupledValue("strainyy"))

{
}

Real
CompXi::computeValue()
{
  Real xi = (_strainxx[_qp] + _strainyy[_qp]) / sqrt( _strainxx[_qp] * _strainxx[_qp] + _strainyy[_qp] * _strainyy[_qp] + 2 * _strainxy[_qp] * _strainxy[_qp] );
  return xi;
}