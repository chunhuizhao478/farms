//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PoroCZMInterfaceKernelPressure.h"

registerMooseObject("farmsApp", PoroCZMInterfaceKernelPressure);

InputParameters
PoroCZMInterfaceKernelPressure::validParams()
{
  InputParameters params = PoroCZMInterfaceKernelPressurebase::validParams();

  params.addClassDescription(
      "CZM Interface kernel to use when using the Small Strain kinematic formulation.");
  return params;
}

PoroCZMInterfaceKernelPressure::PoroCZMInterfaceKernelPressure(const InputParameters & parameters)
  : PoroCZMInterfaceKernelPressurebase(parameters)
{
}


Real
PoroCZMInterfaceKernelPressure::computeDResidualDPressure( const Moose::DGJacobianType & type) const
{
  Real jac;
  if (_component == 0)
    jac = _dtraction_dpressure_global[_qp](_component) - 1;
  else
    jac = _dtraction_dpressure_global[_qp](_component);
    
  switch (type)
  {
    case Moose::ElementElement: // Residual_sign -1  dpf_dpressure sign 1;
      jac *= - _test[_i][_qp] * _phi[_j][_qp];
      break;
    case Moose::ElementNeighbor: // Residual_sign -1   dpf_dpressure sign 1;
      jac *= -_test[_i][_qp] * _phi_neighbor[_j][_qp];
      break;
    case Moose::NeighborElement: // Residual_sign 1   dpf_dpressure sign 1;
      jac *= _test_neighbor[_i][_qp] * _phi[_j][_qp];
      break;
    case Moose::NeighborNeighbor: // Residual_sign 1   dpf_dpressure sign 1;
      jac *= _test_neighbor[_i][_qp] * _phi_neighbor[_j][_qp];
      break;
  }
  return jac;
}

