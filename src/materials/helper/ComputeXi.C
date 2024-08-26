//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ComputeXi.h"

/**
 *  Created by Chunhui Zhao, Jul 17th, 2024
 *  Material used to compute xi of elastic strain
 */
registerMooseObject("farmsApp", ComputeXi);

InputParameters
ComputeXi::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Material used to compute xi of elastic strain");
  return params;
}

ComputeXi::ComputeXi(const InputParameters & parameters)
  : Material(parameters),
    _xi(declareProperty<Real>("xi_initial")),
    _I1(declareProperty<Real>("I1_initial")),
    _I2(declareProperty<Real>("I2_initial")),
    _mechanical_strain(getMaterialProperty<RankTwoTensor>("mechanical_strain"))
{
}

void
ComputeXi::computeQpProperties()
{
    _I1[_qp] = (_mechanical_strain[_qp](0,0)+_mechanical_strain[_qp](1,1)+_mechanical_strain[_qp](2,2));
    _I2[_qp] = _mechanical_strain[_qp](0,0)*_mechanical_strain[_qp](0,0)+_mechanical_strain[_qp](1,1)*_mechanical_strain[_qp](1,1)+_mechanical_strain[_qp](2,2)*_mechanical_strain[_qp](2,2)+2*_mechanical_strain[_qp](0,1)*_mechanical_strain[_qp](0,1)+2*_mechanical_strain[_qp](0,2)*_mechanical_strain[_qp](0,2)+2*_mechanical_strain[_qp](1,2)*_mechanical_strain[_qp](1,2);
    _xi[_qp] = (_mechanical_strain[_qp](0,0)+_mechanical_strain[_qp](1,1)+_mechanical_strain[_qp](2,2))/std::sqrt(_mechanical_strain[_qp](0,0)*_mechanical_strain[_qp](0,0)+_mechanical_strain[_qp](1,1)*_mechanical_strain[_qp](1,1)+_mechanical_strain[_qp](2,2)*_mechanical_strain[_qp](2,2)+2*_mechanical_strain[_qp](0,1)*_mechanical_strain[_qp](0,1)+2*_mechanical_strain[_qp](0,2)*_mechanical_strain[_qp](0,2)+2*_mechanical_strain[_qp](1,2)*_mechanical_strain[_qp](1,2));
}