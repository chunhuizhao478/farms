#include "FarmsInterfaceKernelGlobalz.h"

registerMooseObject("farmsApp", FarmsInterfaceKernelGlobalz);

InputParameters
FarmsInterfaceKernelGlobalz::validParams()
{
  InputParameters params = InterfaceKernel::validParams();
  params.addClassDescription("SlipWeakening Frictional Law Interface Kernel x dir.");
  return params;
}

FarmsInterfaceKernelGlobalz::FarmsInterfaceKernelGlobalz(const InputParameters & parameters)
  : InterfaceKernel(parameters),
    _displacements_plus_global(getMaterialPropertyByName<RealVectorValue>("displacements_plus_global")),
    _displacements_minus_global(getMaterialPropertyByName<RealVectorValue>("displacements_minus_global"))
{
}

Real
FarmsInterfaceKernelGlobalz::computeQpResidual(Moose::DGResidualType type)
{
  Real r = 0;
  switch (type)
  {
    // displacement applied on the primary surface
    case Moose::Element:
      r = ( _displacements_plus_global[_qp](2) ) * _test[_i][_qp];
      break;

    // displacement applied on the secondary surface
    case Moose::Neighbor:
      r = ( _displacements_minus_global[_qp](2) ) * _test_neighbor[_i][_qp];
      break;
  }
  return r;
}

// Real
// FarmsInterfaceKernelGlobalz::computeQpJacobian(Moose::DGJacobianType type)
// {
//   Real jac = 0;
//   return jac;
// }
