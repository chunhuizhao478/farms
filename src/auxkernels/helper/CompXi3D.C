/*
AuxKernel of Passing Variable
*/

#include "CompXi3D.h"

registerMooseObject("farmsApp", CompXi3D);

InputParameters
CompXi3D::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addRequiredCoupledVar("strainxx","strain xx");
  params.addRequiredCoupledVar("strainxy","strain xy");
  params.addRequiredCoupledVar("strainyy","strain yy");
  params.addRequiredCoupledVar("strainxz","strain xz");
  params.addRequiredCoupledVar("strainyz","strain yz");
  params.addRequiredCoupledVar("strainzz","strain zz");
  return params;
}

CompXi3D::CompXi3D(const InputParameters & parameters)
  : AuxKernel(parameters),
  _strainxx(coupledValue("strainxx")),
  _strainxy(coupledValue("strainxy")),
  _strainyy(coupledValue("strainyy")),
  _strainxz(coupledValue("strainxz")),
  _strainyz(coupledValue("strainyz")),
  _strainzz(coupledValue("strainzz"))
{
}

Real
CompXi3D::computeValue()
{
  Real xi = (_strainxx[_qp] + _strainyy[_qp] + _strainzz[_qp]) / sqrt( _strainxx[_qp] * _strainxx[_qp] + _strainyy[_qp] * _strainyy[_qp] + _strainzz[_qp] * _strainzz[_qp] + 2 * _strainxy[_qp] * _strainxy[_qp] + 2 * _strainxz[_qp] * _strainxz[_qp] + 2 * _strainyz[_qp] * _strainyz[_qp] );
  return xi;
}