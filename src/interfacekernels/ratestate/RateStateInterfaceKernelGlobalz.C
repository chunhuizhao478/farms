#include "RateStateInterfaceKernelGlobalz.h"

registerMooseObject("farmsApp", RateStateInterfaceKernelGlobalz);

InputParameters
RateStateInterfaceKernelGlobalz::validParams()
{
  InputParameters params = InterfaceKernel::validParams();
  params.addClassDescription("Rate-and-State Frictional Law Interface Kernel x dir.");
  return params;
}

RateStateInterfaceKernelGlobalz::RateStateInterfaceKernelGlobalz(const InputParameters & parameters)
  : InterfaceKernel(parameters),
  _disp_dip_plus(getMaterialPropertyByName<Real>("alongfaultdisp_dip_plus")),
  _disp_dip_minus(getMaterialPropertyByName<Real>("alongfaultdisp_dip_minus"))
{
}

Real
RateStateInterfaceKernelGlobalz::computeQpResidual(Moose::DGResidualType type)
{

  //!rot is NOT implemented!

  Real r = 0;
  switch (type)
  {
    // displacement applied on the primary surface
    case Moose::Element:
      r = ( _disp_dip_plus[_qp] ) * _test[_i][_qp];
      break;

    // displacement applied on the secondary surface
    case Moose::Neighbor:
      r = ( _disp_dip_minus[_qp] ) * _test_neighbor[_i][_qp];
      break;
  }
  return r;
}

// Real
// RateStateInterfaceKernelGlobalz::computeQpJacobian(Moose::DGJacobianType type)
// {
//   Real jac = 0;
//   return jac;
// }
