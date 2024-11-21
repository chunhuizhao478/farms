#include "RateStateInterfaceKernelGlobaly.h"

registerMooseObject("farmsApp", RateStateInterfaceKernelGlobaly);

InputParameters
RateStateInterfaceKernelGlobaly::validParams()
{
  InputParameters params = InterfaceKernel::validParams();
  params.addClassDescription("Rate-and-State Frictional Law Interface Kernel y dir.");
  params.addRequiredCoupledVar("x_var", "x-displacement variable");
  return params;
}

RateStateInterfaceKernelGlobaly::RateStateInterfaceKernelGlobaly(const InputParameters & parameters)
  : InterfaceKernel(parameters),
    _disp_normal_plus(getMaterialPropertyByName<Real>("alongfaultdisp_normal_plus")),
    _disp_normal_minus(getMaterialPropertyByName<Real>("alongfaultdisp_normal_minus")),
    _dalongfaultdisp_normal_plus_tplusdt_du_strike_plus(getMaterialPropertyByName<Real>("dalongfaultdisp_normal_plus_tplusdt_du_strike_plus")),
    _dalongfaultdisp_normal_plus_tplusdt_du_strike_minus(getMaterialPropertyByName<Real>("dalongfaultdisp_normal_plus_tplusdt_du_strike_minus")),
    _dalongfaultdisp_normal_plus_tplusdt_du_normal_plus(getMaterialPropertyByName<Real>("dalongfaultdisp_normal_plus_tplusdt_du_normal_plus")),
    _dalongfaultdisp_normal_plus_tplusdt_du_normal_minus(getMaterialPropertyByName<Real>("dalongfaultdisp_normal_plus_tplusdt_du_normal_minus")),
    _dalongfaultdisp_normal_minus_tplusdt_du_strike_plus(getMaterialPropertyByName<Real>("dalongfaultdisp_normal_minus_tplusdt_du_strike_plus")),
    _dalongfaultdisp_normal_minus_tplusdt_du_strike_minus(getMaterialPropertyByName<Real>("dalongfaultdisp_normal_minus_tplusdt_du_strike_minus")),
    _dalongfaultdisp_normal_minus_tplusdt_du_normal_plus(getMaterialPropertyByName<Real>("dalongfaultdisp_normal_minus_tplusdt_du_normal_plus")),
    _dalongfaultdisp_normal_minus_tplusdt_du_normal_minus(getMaterialPropertyByName<Real>("dalongfaultdisp_normal_minus_tplusdt_du_normal_minus")),
    _x_var(coupled("x_var"))
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

Real
RateStateInterfaceKernelGlobaly::computeQpJacobian(Moose::DGJacobianType type)
{
  Real jac = 0;
  switch (type)
  {
    case Moose::ElementElement:      // d(Element y-Residual)/d(Element y-DOF)
      jac = -1.0 * _dalongfaultdisp_normal_plus_tplusdt_du_normal_plus[_qp]  * _test[_i][_qp] * _phi[_j][_qp];
      break;
      
    case Moose::ElementNeighbor:     // d(Element y-Residual)/d(Neighbor y-DOF)
      jac = -1.0 * _dalongfaultdisp_normal_plus_tplusdt_du_normal_minus[_qp]  * _test[_i][_qp] * _phi_neighbor[_j][_qp];
      break;
      
    case Moose::NeighborElement:     // d(Neighbor y-Residual)/d(Element y-DOF)
      jac = -1.0 * _dalongfaultdisp_normal_minus_tplusdt_du_normal_plus[_qp]  * _test_neighbor[_i][_qp] * _phi[_j][_qp];
      break;
      
    case Moose::NeighborNeighbor:    // d(Neighbor y-Residual)/d(Neighbor y-DOF)
      jac = -1.0 * _dalongfaultdisp_normal_minus_tplusdt_du_normal_minus[_qp]  * _test_neighbor[_i][_qp] * _phi_neighbor[_j][_qp];
      break;
  }
  return jac;
}

Real
RateStateInterfaceKernelGlobaly::computeQpOffDiagJacobian(Moose::DGJacobianType type, unsigned int jvar)
{
  Real jac = 0;
  
  if (jvar == _x_var)  // Coupling with x direction
  {
    switch (type)
    {
      case Moose::ElementElement:      // d(Element y-Residual)/d(Element x-DOF)
        jac = -1.0 * _dalongfaultdisp_normal_plus_tplusdt_du_strike_plus[_qp]  * _test[_i][_qp] * _phi[_j][_qp];
        break;
        
      case Moose::ElementNeighbor:     // d(Element y-Residual)/d(Neighbor x-DOF)
        jac = -1.0 * _dalongfaultdisp_normal_plus_tplusdt_du_strike_minus[_qp]  * _test[_i][_qp] * _phi_neighbor[_j][_qp];
        break;
        
      case Moose::NeighborElement:     // d(Neighbor y-Residual)/d(Element x-DOF)
        jac = -1.0 * _dalongfaultdisp_normal_minus_tplusdt_du_strike_plus[_qp]  * _test_neighbor[_i][_qp] * _phi[_j][_qp];
        break;
        
      case Moose::NeighborNeighbor:    // d(Neighbor y-Residual)/d(Neighbor x-DOF)
        jac = -1.0 * _dalongfaultdisp_normal_minus_tplusdt_du_strike_minus[_qp]  * _test_neighbor[_i][_qp] * _phi_neighbor[_j][_qp];
        break;
    }
  }
  return jac;
}