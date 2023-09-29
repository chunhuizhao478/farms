#include "RateStateInterfaceKernelGlobalxdev.h"

registerMooseObject("farmsApp", RateStateInterfaceKernelGlobalxdev);

InputParameters
RateStateInterfaceKernelGlobalxdev::validParams()
{
  InputParameters params = InterfaceKernel::validParams();
  params.addClassDescription("Rate-and-State Frictional Law Interface Kernel x dir.");
  params.addRequiredCoupledVar("disp_strike_plus","along fault displacement plus side in x dir");
  params.addRequiredCoupledVar("disp_strike_minus","along fault displacement minus side in x dir");
  return params;
}

RateStateInterfaceKernelGlobalxdev::RateStateInterfaceKernelGlobalxdev(const InputParameters & parameters)
  : InterfaceKernel(parameters),
  _disp_strike_plus(coupledValue("disp_strike_plus")),
  _disp_strike_minus(coupledValue("disp_strike_minus"))
{
}

Real
RateStateInterfaceKernelGlobalxdev::computeQpResidual(Moose::DGResidualType type)
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
// RateStateInterfaceKernelGlobalxdev::computeQpJacobian(Moose::DGJacobianType type)
// {
//   Real jac = 0;
//   return jac;
// }
