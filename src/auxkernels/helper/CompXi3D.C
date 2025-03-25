/*
AuxKernel of Passing Variable
*/

#include "CompXi3D.h"

registerMooseObject("farmsApp", CompXi3D);

InputParameters
CompXi3D::validParams()
{
  InputParameters params = AuxKernel::validParams();
  return params;
}

CompXi3D::CompXi3D(const InputParameters & parameters)
  : AuxKernel(parameters),
  _mechanical_strain(getMaterialProperty<RankTwoTensor>("mechanical_strain"))
{
}

Real
CompXi3D::computeValue()
{

  const Real epsilon = 1e-12;
  Real I1 = epsilon + _mechanical_strain[_qp](0,0) + _mechanical_strain[_qp](1,1) + _mechanical_strain[_qp](2,2);
  Real I2 = epsilon + _mechanical_strain[_qp](0,0) * _mechanical_strain[_qp](0,0) + _mechanical_strain[_qp](1,1) * _mechanical_strain[_qp](1,1) + _mechanical_strain[_qp](2,2) * _mechanical_strain[_qp](2,2) + 2 * _mechanical_strain[_qp](1,2) * _mechanical_strain[_qp](1,2) + 2 * _mechanical_strain[_qp](0,1) * _mechanical_strain[_qp](0,1) + 2 * _mechanical_strain[_qp](0,2) * _mechanical_strain[_qp](0,2);
  Real xi = I1 / std::sqrt(I2);

  //Real xi = (_mechanical_strain[_qp](0,0) + _mechanical_strain[_qp](1,1) + _mechanical_strain[_qp](2,2)) / sqrt( _mechanical_strain[_qp](0,0) * _mechanical_strain[_qp](0,0) + _mechanical_strain[_qp](1,1) * _mechanical_strain[_qp](1,1) + _mechanical_strain[_qp](2,2) * _mechanical_strain[_qp](2,2) + 2 * _mechanical_strain[_qp](0,1) * _mechanical_strain[_qp](0,1) + 2 * _mechanical_strain[_qp](0,2) * _mechanical_strain[_qp](0,2) + 2 * _mechanical_strain[_qp](1,2) * _mechanical_strain[_qp](1,2) );
  return xi;
}