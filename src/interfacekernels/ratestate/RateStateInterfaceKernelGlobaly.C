#include "RateStateInterfaceKernelGlobaly.h"

registerMooseObject("farmsApp", RateStateInterfaceKernelGlobaly);

InputParameters
RateStateInterfaceKernelGlobaly::validParams()
{
  InputParameters params = InterfaceKernel::validParams();
  params.addClassDescription("Rate-and-State Frictional Law Interface Kernel y dir.");
  return params;
}

RateStateInterfaceKernelGlobaly::RateStateInterfaceKernelGlobaly(const InputParameters & parameters)
  : InterfaceKernel(parameters),
  _disp_normal_plus(getMaterialPropertyByName<Real>("alongfaultdisp_normal_plus")),
  _disp_normal_minus(getMaterialPropertyByName<Real>("alongfaultdisp_normal_minus"))
{
}

Real
RateStateInterfaceKernelGlobaly::computeQpResidual(Moose::DGResidualType type)
{

  //!rot is NOT implemented!

  Real r = 0;
  switch (type)
  {
    // displacement applied on the primary surface
    case Moose::Element:
      r = ( -1 * _disp_normal_plus[_qp] ) * _test[_i][_qp];
      break;

    // displacement applied on the secondary surface
    case Moose::Neighbor:
      r = ( -1 * _disp_normal_minus[_qp] ) * _test_neighbor[_i][_qp];
      break;
  }
  return r;
}

// Real
// RateStateInterfaceKernelGlobaly::computeQpJacobian(Moose::DGJacobianType type)
// {
//   Real jac = 0;
//   return jac;
// }
