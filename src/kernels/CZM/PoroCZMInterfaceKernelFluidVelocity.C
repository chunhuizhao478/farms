//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PoroCZMInterfaceKernelFluidVelocity.h"

registerMooseObject("farmsApp", PoroCZMInterfaceKernelFluidVelocity);

InputParameters
PoroCZMInterfaceKernelFluidVelocity::validParams()
{
  InputParameters params = PoroCZMInterfaceKernelFluxbase::validParams();

  params.addClassDescription(
      "CZM Interface kernel to use when using the Small Strain kinematic formulation.");
  return params;
}

PoroCZMInterfaceKernelFluidVelocity::PoroCZMInterfaceKernelFluidVelocity(const InputParameters & parameters)
  : PoroCZMInterfaceKernelFluxbase(parameters)
{
}


Real
PoroCZMInterfaceKernelFluidVelocity::computeDResidualDFlux(
    const unsigned int & component_j, const Moose::DGJacobianType & type) const
{
  Real jac;
  if (_component == 0)
    jac = _dtraction_djump_global_vf[_qp](_component, component_j) - _dpressure_djump_global_vf[_qp](component_j);
  else
    jac = _dtraction_djump_global_vf[_qp](_component, component_j);
    
  switch (type)
  {
    case Moose::ElementElement: // Residual_sign -1  ddeltaV_dfluidvel sign -1;
      jac *= _test[_i][_qp] * _fluid_vars[component_j]->phiFace()[_j][_qp];
      break;
    case Moose::ElementNeighbor: // Residual_sign -1  ddeltaV_dfluidvel sign 1;
      jac *= -_test[_i][_qp] * _fluid_vars[component_j]->phiFaceNeighbor()[_j][_qp];
      break;
    case Moose::NeighborElement: // Residual_sign 1  ddeltaV_dfluidvel sign -1;
      jac *= -_test_neighbor[_i][_qp] * _fluid_vars[component_j]->phiFace()[_j][_qp];
      break;
    case Moose::NeighborNeighbor: // Residual_sign 1  ddeltaV_dfluidvel sign 1;
      jac *= _test_neighbor[_i][_qp] * _fluid_vars[component_j]->phiFaceNeighbor()[_j][_qp];
      break;
  }
  return jac;
}
