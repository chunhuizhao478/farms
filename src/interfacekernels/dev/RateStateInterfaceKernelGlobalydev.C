#include "RateStateInterfaceKernelGlobalydev.h"

registerMooseObject("farmsApp", RateStateInterfaceKernelGlobalydev);

InputParameters
RateStateInterfaceKernelGlobalydev::validParams()
{
  InputParameters params = InterfaceKernel::validParams();
  params.addClassDescription("Rate-and-State Frictional Law Interface Kernel y dir.");
  params.addRequiredCoupledVar("disp_normal_plus","along fault displacement plus side in y dir");
  params.addRequiredCoupledVar("disp_normal_minus","along fault displacement minus side in y dir");
  return params;
}

RateStateInterfaceKernelGlobalydev::RateStateInterfaceKernelGlobalydev(const InputParameters & parameters)
  : InterfaceKernel(parameters),
  _disp_normal_plus(coupledValue("disp_normal_plus")),
  _disp_normal_minus(coupledValue("disp_normal_minus"))
{
}

Real
RateStateInterfaceKernelGlobalydev::computeQpResidual(Moose::DGResidualType type)
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
// RateStateInterfaceKernelGlobalydev::computeQpJacobian(Moose::DGJacobianType type)
// {
//   Real jac = 0;
//   return jac;
// }
