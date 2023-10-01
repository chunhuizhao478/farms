#include "RateStateInterfaceKernelGlobalx.h"

registerMooseObject("farmsApp", RateStateInterfaceKernelGlobalx);

InputParameters
RateStateInterfaceKernelGlobalx::validParams()
{
  InputParameters params = InterfaceKernel::validParams();
  params.addClassDescription("Rate-and-State Frictional Law Interface Kernel x dir.");
  return params;
}

RateStateInterfaceKernelGlobalx::RateStateInterfaceKernelGlobalx(const InputParameters & parameters)
  : InterfaceKernel(parameters),
  _disp_strike_plus(getMaterialPropertyByName<Real>("alongfaultdisp_strike_plus")),
  _disp_strike_minus(getMaterialPropertyByName<Real>("alongfaultdisp_strike_minus"))
{
}

Real
RateStateInterfaceKernelGlobalx::computeQpResidual(Moose::DGResidualType type)
{

  //!rot is NOT implemented!

  Real r = 0;
  switch (type)
  {
    // displacement applied on the primary surface
    case Moose::Element:
      r = ( _disp_strike_plus[_qp] ) * _test[_i][_qp];
      break;

    // displacement applied on the secondary surface
    case Moose::Neighbor:
      r = ( _disp_strike_minus[_qp] ) * _test_neighbor[_i][_qp];
      break;
  }
  return r;
}

// Real
// RateStateInterfaceKernelGlobalx::computeQpJacobian(Moose::DGJacobianType type)
// {
//   Real jac = 0;
//   return jac;
// }
