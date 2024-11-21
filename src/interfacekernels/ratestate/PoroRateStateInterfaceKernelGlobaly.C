#include "PoroRateStateInterfaceKernelGlobaly.h"

registerMooseObject("farmsApp", PoroRateStateInterfaceKernelGlobaly);

InputParameters
PoroRateStateInterfaceKernelGlobaly::validParams()
{
  InputParameters params = InterfaceKernel::validParams();
  params.addClassDescription("Rate-and-State Frictional Law Interface Kernel y dir.");
  params.addRequiredCoupledVar("x_var", "x-displacement variable");
  params.addRequiredCoupledVar("vfx_var", "x-fluid velocity variable");
  params.addRequiredCoupledVar("vfy_var", "y-fluid velocity variable");
  params.addRequiredCoupledVar("pressure_var", "pore pressure variable");
  return params;
}

PoroRateStateInterfaceKernelGlobaly::PoroRateStateInterfaceKernelGlobaly(const InputParameters & parameters)
  : InterfaceKernel(parameters),
    _base_name(isParamValid("base_name") && !getParam<std::string>("base_name").empty()
               ? getParam<std::string>("base_name") + "_"
               : ""),
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
    _dalongfaultdisp_normal_plus_tplusdt_dvf_strike_plus(getMaterialPropertyByName<Real>("dalongfaultdisp_normal_plus_tplusdt_dvf_strike_plus")),
    _dalongfaultdisp_normal_plus_tplusdt_dvf_strike_minus(getMaterialPropertyByName<Real>("dalongfaultdisp_normal_plus_tplusdt_dvf_strike_minus")),
    _dalongfaultdisp_normal_plus_tplusdt_dvf_normal_plus(getMaterialPropertyByName<Real>("dalongfaultdisp_normal_plus_tplusdt_dvf_normal_plus")),
    _dalongfaultdisp_normal_plus_tplusdt_dvf_normal_minus(getMaterialPropertyByName<Real>("dalongfaultdisp_normal_plus_tplusdt_dvf_normal_minus")),
    _dalongfaultdisp_normal_minus_tplusdt_dvf_strike_plus(getMaterialPropertyByName<Real>("dalongfaultdisp_normal_minus_tplusdt_dvf_strike_plus")),
    _dalongfaultdisp_normal_minus_tplusdt_dvf_strike_minus(getMaterialPropertyByName<Real>("dalongfaultdisp_normal_minus_tplusdt_dvf_strike_minus")),
    _dalongfaultdisp_normal_minus_tplusdt_dvf_normal_plus(getMaterialPropertyByName<Real>("dalongfaultdisp_normal_minus_tplusdt_dvf_normal_plus")),
    _dalongfaultdisp_normal_minus_tplusdt_dvf_normal_minus(getMaterialPropertyByName<Real>("dalongfaultdisp_normal_minus_tplusdt_dvf_normal_minus")),
    _dalongfaultdisp_normal_plus_tplusdt_dp_plus(getMaterialPropertyByName<Real>("dalongfaultdisp_normal_plus_tplusdt_dp_plus")),
    _dalongfaultdisp_normal_plus_tplusdt_dp_minus(getMaterialPropertyByName<Real>("dalongfaultdisp_normal_plus_tplusdt_dp_minus")),
    _dalongfaultdisp_normal_minus_tplusdt_dp_plus(getMaterialPropertyByName<Real>("dalongfaultdisp_normal_minus_tplusdt_dp_plus")),
    _dalongfaultdisp_normal_minus_tplusdt_dp_minus(getMaterialPropertyByName<Real>("dalongfaultdisp_normal_minus_tplusdt_dp_minus")),
    _x_var(coupled("x_var")),
    _vfx_var(coupled("vfx_var")),
    _vfy_var(coupled("vfy_var")),
    _pressure_var(coupled("pressure_var")),
    _vfx_var_data(*getVar("vfx_var", 0)),
    _vfy_var_data(*getVar("vfy_var", 0)),
    _pressure_var_data(*getVar("pressure_var", 0)),
    _phi_vfx(_assembly.phiFace(_vfx_var_data)),
    _phi_vfx_neighbor(_assembly.phiFaceNeighbor(_vfx_var_data)),
    _phi_vfy(_assembly.phiFace(_vfy_var_data)),
    _phi_vfy_neighbor(_assembly.phiFaceNeighbor(_vfy_var_data)),
    _phi_pressure(_assembly.phiFace(_pressure_var_data)),
    _phi_pressure_neighbor(_assembly.phiFaceNeighbor(_pressure_var_data))
{
}

Real
PoroRateStateInterfaceKernelGlobaly::computeQpResidual(Moose::DGResidualType type)
{
  Real r = 0;
  switch (type)
  {
    case Moose::Element:
      r = (-1 * _disp_normal_plus[_qp]) * _test[_i][_qp];
      break;

    case Moose::Neighbor:
      r = (-1 * _disp_normal_minus[_qp]) * _test_neighbor[_i][_qp];
      break;
  }
  return r;
}

Real
PoroRateStateInterfaceKernelGlobaly::computeQpJacobian(Moose::DGJacobianType type)
{
  Real jac = 0;
  switch (type)
  {
    case Moose::ElementElement:
      jac = -1.0 * _dalongfaultdisp_normal_plus_tplusdt_du_normal_plus[_qp] * _test[_i][_qp] * _phi[_j][_qp];
      break;
      
    case Moose::ElementNeighbor:
      jac = -1.0 * _dalongfaultdisp_normal_plus_tplusdt_du_normal_minus[_qp] * _test[_i][_qp] * _phi_neighbor[_j][_qp];
      break;
      
    case Moose::NeighborElement:
      jac = -1.0 * _dalongfaultdisp_normal_minus_tplusdt_du_normal_plus[_qp] * _test_neighbor[_i][_qp] * _phi[_j][_qp];
      break;
      
    case Moose::NeighborNeighbor:
      jac = -1.0 * _dalongfaultdisp_normal_minus_tplusdt_du_normal_minus[_qp] * _test_neighbor[_i][_qp] * _phi_neighbor[_j][_qp];
      break;
  }
  return jac;
}

Real
PoroRateStateInterfaceKernelGlobaly::computeQpOffDiagJacobian(Moose::DGJacobianType type, unsigned int jvar)
{
  Real jac = 0;
  
  if (jvar == _x_var)  // Coupling with x direction (second order)
  {
    switch (type)
    {
      case Moose::ElementElement:
        jac = -1.0 * _dalongfaultdisp_normal_plus_tplusdt_du_strike_plus[_qp] * _test[_i][_qp] * _phi[_j][_qp];
        break;
        
      case Moose::ElementNeighbor:
        jac = -1.0 * _dalongfaultdisp_normal_plus_tplusdt_du_strike_minus[_qp] * _test[_i][_qp] * _phi_neighbor[_j][_qp];
        break;
        
      case Moose::NeighborElement:
        jac = -1.0 * _dalongfaultdisp_normal_minus_tplusdt_du_strike_plus[_qp] * _test_neighbor[_i][_qp] * _phi[_j][_qp];
        break;
        
      case Moose::NeighborNeighbor:
        jac = -1.0 * _dalongfaultdisp_normal_minus_tplusdt_du_strike_minus[_qp] * _test_neighbor[_i][_qp] * _phi_neighbor[_j][_qp];
        break;
    }
  }
  else if (jvar == _vfx_var) // Coupling with x fluid velocity (first order)
  {
    switch (type)
    {
      case Moose::ElementElement:
        jac = -1.0 * _dalongfaultdisp_normal_plus_tplusdt_dvf_strike_plus[_qp] * _test[_i][_qp] * _phi_vfx[_j][_qp];
        break;
        
      case Moose::ElementNeighbor:
        jac = -1.0 * _dalongfaultdisp_normal_plus_tplusdt_dvf_strike_minus[_qp] * _test[_i][_qp] * _phi_vfx_neighbor[_j][_qp];
        break;
        
      case Moose::NeighborElement:
        jac = -1.0 * _dalongfaultdisp_normal_minus_tplusdt_dvf_strike_plus[_qp] * _test_neighbor[_i][_qp] * _phi_vfx[_j][_qp];
        break;
        
      case Moose::NeighborNeighbor:
        jac = -1.0 * _dalongfaultdisp_normal_minus_tplusdt_dvf_strike_minus[_qp] * _test_neighbor[_i][_qp] * _phi_vfx_neighbor[_j][_qp];
        break;
    }
  }
  else if (jvar == _vfy_var) // Coupling with y fluid velocity (first order)
  {
    switch (type)
    {
      case Moose::ElementElement:
        jac = -1.0 * _dalongfaultdisp_normal_plus_tplusdt_dvf_normal_plus[_qp] * _test[_i][_qp] * _phi_vfy[_j][_qp];
        break;
        
      case Moose::ElementNeighbor:
        jac = -1.0 * _dalongfaultdisp_normal_plus_tplusdt_dvf_normal_minus[_qp] * _test[_i][_qp] * _phi_vfy_neighbor[_j][_qp];
        break;
        
      case Moose::NeighborElement:
        jac = -1.0 * _dalongfaultdisp_normal_minus_tplusdt_dvf_normal_plus[_qp] * _test_neighbor[_i][_qp] * _phi_vfy[_j][_qp];
        break;
        
      case Moose::NeighborNeighbor:
        jac = -1.0 * _dalongfaultdisp_normal_minus_tplusdt_dvf_normal_minus[_qp] * _test_neighbor[_i][_qp] * _phi_vfy_neighbor[_j][_qp];
        break;
    }
  }
  else if (jvar == _pressure_var) // Coupling with pressure (first order)
  {
    switch (type)
    {
      case Moose::ElementElement:
        jac = -1.0 * _dalongfaultdisp_normal_plus_tplusdt_dp_plus[_qp] * _test[_i][_qp] * _phi_pressure[_j][_qp];
        break;
        
      case Moose::ElementNeighbor:
        jac = -1.0 * _dalongfaultdisp_normal_plus_tplusdt_dp_minus[_qp] * _test[_i][_qp] * _phi_pressure_neighbor[_j][_qp];
        break;
        
      case Moose::NeighborElement:
        jac = -1.0 * _dalongfaultdisp_normal_minus_tplusdt_dp_plus[_qp] * _test_neighbor[_i][_qp] * _phi_pressure[_j][_qp];
        break;
        
      case Moose::NeighborNeighbor:
        jac = -1.0 * _dalongfaultdisp_normal_minus_tplusdt_dp_minus[_qp] * _test_neighbor[_i][_qp] * _phi_pressure_neighbor[_j][_qp];
        break;
    }
  }
  return jac;
}