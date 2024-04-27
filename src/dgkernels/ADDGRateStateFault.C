//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADDGRateStateFault.h"

// MOOSE includes
#include "MooseVariableFE.h"

#include "libmesh/utility.h"

registerMooseObject("farmsApp", ADDGRateStateFault);

InputParameters
ADDGRateStateFault::validParams()
{
  InputParameters params = ADDGKernel::validParams();
  // See header file for sigma and epsilon
  params.addRequiredParam<Real>("sigma", "sigma");
  params.addRequiredParam<Real>("epsilon", "epsilon");
  params.addParam<MaterialPropertyName>(
      "diff", 1., "The diffusion (or thermal conductivity or viscosity) coefficient.");
  params.addClassDescription("DG kernel for diffusion operator");
  return params;
}

ADDGRateStateFault::ADDGRateStateFault(const InputParameters & parameters)
  : ADDGKernel(parameters),
    _epsilon(getParam<Real>("epsilon")),
    _sigma(getParam<Real>("sigma")),
    _diff(getADMaterialProperty<Real>("diff")),
    _diff_neighbor(getNeighborADMaterialProperty<Real>("diff")),
    _elasticity(getMaterialProperty<Real>("elasticity")),
    _slip(getADMaterialProperty<Real>("slip"))
{
}

ADReal
ADDGRateStateFault::computeQpResidual(Moose::DGResidualType type)
{
  ADReal r = 0.0;

  const int elem_b_order = std::max(libMesh::Order(1), _var.order());
  const Real h_elem =
      _current_elem_volume / _current_side_volume * 1.0 / Utility::pow<2>(elem_b_order);

  switch (type)
  {
    case Moose::Element:
      r -= _sigma / h_elem * _test[_i][_qp] * _slip[_qp];
      r += 0.5 * _slip[_qp] * _elasticity[_qp] * _grad_test[_i][_qp] * _normals[_qp];
      break;

    case Moose::Neighbor:
      r += _sigma / h_elem * _test_neighbor[_i][_qp] * _slip[_qp];
      r += 0.5 * _slip[_qp] * _elasticity[_qp] * _grad_test_neighbor[_i][_qp] * _normals[_qp];
      break;
  }

  return r;
}