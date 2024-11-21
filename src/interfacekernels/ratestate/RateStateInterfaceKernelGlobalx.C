#include "RateStateInterfaceKernelGlobalx.h"

registerMooseObject("farmsApp", RateStateInterfaceKernelGlobalx);

InputParameters
RateStateInterfaceKernelGlobalx::validParams()
{
  InputParameters params = InterfaceKernel::validParams();
  params.addClassDescription("Rate-and-State Frictional Law Interface Kernel x dir.");
  params.addRequiredCoupledVar("y_var", "y-displacement variable");
  return params;
}

RateStateInterfaceKernelGlobalx::RateStateInterfaceKernelGlobalx(const InputParameters & parameters)
  : InterfaceKernel(parameters),
    _disp_strike_plus(getMaterialPropertyByName<Real>("alongfaultdisp_strike_plus")),
    _disp_strike_minus(getMaterialPropertyByName<Real>("alongfaultdisp_strike_minus")),
    _dalongfaultdisp_strike_plus_tplusdt_du_strike_plus(getMaterialPropertyByName<Real>("dalongfaultdisp_strike_plus_tplusdt_du_strike_plus")),
    _dalongfaultdisp_strike_plus_tplusdt_du_strike_minus(getMaterialPropertyByName<Real>("dalongfaultdisp_strike_plus_tplusdt_du_strike_minus")),
    _dalongfaultdisp_strike_plus_tplusdt_du_normal_plus(getMaterialPropertyByName<Real>("dalongfaultdisp_strike_plus_tplusdt_du_normal_plus")),
    _dalongfaultdisp_strike_plus_tplusdt_du_normal_minus(getMaterialPropertyByName<Real>("dalongfaultdisp_strike_plus_tplusdt_du_normal_minus")),
    _dalongfaultdisp_strike_minus_tplusdt_du_strike_plus(getMaterialPropertyByName<Real>("dalongfaultdisp_strike_minus_tplusdt_du_strike_plus")),
    _dalongfaultdisp_strike_minus_tplusdt_du_strike_minus(getMaterialPropertyByName<Real>("dalongfaultdisp_strike_minus_tplusdt_du_strike_minus")),
    _dalongfaultdisp_strike_minus_tplusdt_du_normal_plus(getMaterialPropertyByName<Real>("dalongfaultdisp_strike_minus_tplusdt_du_normal_plus")),
    _dalongfaultdisp_strike_minus_tplusdt_du_normal_minus(getMaterialPropertyByName<Real>("dalongfaultdisp_strike_minus_tplusdt_du_normal_minus")),
    _y_var(coupled("y_var"))
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

Real
RateStateInterfaceKernelGlobalx::computeQpJacobian(Moose::DGJacobianType type)
{
  Real jac = 0;
  switch (type)
  {
    case Moose::ElementElement:      // d(Element x-Residual)/d(Element x-DOF)
      jac = _dalongfaultdisp_strike_plus_tplusdt_du_strike_plus[_qp]  * _test[_i][_qp] * _phi[_j][_qp];
      break;
      
    case Moose::ElementNeighbor:     // d(Element x-Residual)/d(Neighbor x-DOF)
      jac = _dalongfaultdisp_strike_plus_tplusdt_du_strike_minus[_qp]  * _test[_i][_qp] * _phi_neighbor[_j][_qp];
      break;
      
    case Moose::NeighborElement:     // d(Neighbor x-Residual)/d(Element x-DOF)
      jac = _dalongfaultdisp_strike_minus_tplusdt_du_strike_plus[_qp]  * _test_neighbor[_i][_qp] * _phi[_j][_qp];
      break;
      
    case Moose::NeighborNeighbor:    // d(Neighbor x-Residual)/d(Neighbor x-DOF)
      jac = _dalongfaultdisp_strike_minus_tplusdt_du_strike_minus[_qp]  * _test_neighbor[_i][_qp] * _phi_neighbor[_j][_qp];
      break;
  }
  return jac;
}

Real
RateStateInterfaceKernelGlobalx::computeQpOffDiagJacobian(Moose::DGJacobianType type, unsigned int jvar)
{
  Real jac = 0;
  
  if (jvar == _y_var)  // Coupling with y direction
  {
    switch (type)
    {
      case Moose::ElementElement:      // d(Element x-Residual)/d(Element y-DOF)
        jac = _dalongfaultdisp_strike_plus_tplusdt_du_normal_plus[_qp]  * _test[_i][_qp] * _phi[_j][_qp];
        break;
        
      case Moose::ElementNeighbor:     // d(Element x-Residual)/d(Neighbor y-DOF)
        jac = _dalongfaultdisp_strike_plus_tplusdt_du_normal_minus[_qp]  * _test[_i][_qp] * _phi_neighbor[_j][_qp];
        break;
        
      case Moose::NeighborElement:     // d(Neighbor x-Residual)/d(Element y-DOF)
        jac = _dalongfaultdisp_strike_minus_tplusdt_du_normal_plus[_qp]  * _test_neighbor[_i][_qp] * _phi[_j][_qp];
        break;
        
      case Moose::NeighborNeighbor:    // d(Neighbor x-Residual)/d(Neighbor y-DOF)
        jac = _dalongfaultdisp_strike_minus_tplusdt_du_normal_minus[_qp]  * _test_neighbor[_i][_qp] * _phi_neighbor[_j][_qp];
        break;
    }
  }
  return jac;
}
