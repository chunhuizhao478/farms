//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADComputeXi.h"

/**
 *  Created by Chunhui Zhao, Jul 17th, 2024
 *  ADMaterial used to compute xi of elastic strain
 */
registerMooseObject("farmsApp", ADComputeXi);

InputParameters
ADComputeXi::validParams()
{
  InputParameters params = ADMaterial::validParams();
  params.addClassDescription("ADMaterial used to compute xi of elastic strain");
  return params;
}

ADComputeXi::ADComputeXi(const InputParameters & parameters)
  : ADMaterial(parameters),
    _xi(declareADProperty<Real>("xi_initial")),
    _mechanical_strain(getADMaterialProperty<RankTwoTensor>("mechanical_strain"))
{
}

void
ADComputeXi::computeQpProperties()
{
    _xi[_qp] = (_mechanical_strain[_qp](0,0)+_mechanical_strain[_qp](1,1)+_mechanical_strain[_qp](2,2))/std::sqrt(_mechanical_strain[_qp](0,0)*_mechanical_strain[_qp](0,0)+_mechanical_strain[_qp](1,1)*_mechanical_strain[_qp](1,1)+_mechanical_strain[_qp](2,2)*_mechanical_strain[_qp](2,2)+2*_mechanical_strain[_qp](0,1)*_mechanical_strain[_qp](0,1)+2*_mechanical_strain[_qp](0,2)*_mechanical_strain[_qp](0,2)+2*_mechanical_strain[_qp](1,2)*_mechanical_strain[_qp](1,2));
}