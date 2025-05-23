#include "PoroRateStateInterfaceKernelGlobalx.h"

registerMooseObject("farmsApp", PoroRateStateInterfaceKernelGlobalx);

InputParameters
PoroRateStateInterfaceKernelGlobalx::validParams()
{
  InputParameters params = InterfaceKernel::validParams();
  params.addClassDescription("Rate-and-State Frictional Law Interface Kernel x dir.");
  params.addRequiredCoupledVar("y_var", "y-displacement variable");
  params.addRequiredCoupledVar("vfx_var", "x-fluid velocity variable");
  params.addRequiredCoupledVar("vfy_var", "y-fluid velocity variable");
  params.addRequiredCoupledVar("pressure_var", "pore pressure variable");
  return params;
}

PoroRateStateInterfaceKernelGlobalx::PoroRateStateInterfaceKernelGlobalx(const InputParameters & parameters)
  : InterfaceKernel(parameters),
    _base_name(isParamValid("base_name") && !getParam<std::string>("base_name").empty()
               ? getParam<std::string>("base_name") + "_"
               : ""),
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
    _dalongfaultdisp_strike_plus_tplusdt_dvf_strike_plus(getMaterialPropertyByName<Real>("dalongfaultdisp_strike_plus_tplusdt_dvf_strike_plus")),
    _dalongfaultdisp_strike_plus_tplusdt_dvf_strike_minus(getMaterialPropertyByName<Real>("dalongfaultdisp_strike_plus_tplusdt_dvf_strike_minus")),
    _dalongfaultdisp_strike_plus_tplusdt_dvf_normal_plus(getMaterialPropertyByName<Real>("dalongfaultdisp_strike_plus_tplusdt_dvf_normal_plus")),
    _dalongfaultdisp_strike_plus_tplusdt_dvf_normal_minus(getMaterialPropertyByName<Real>("dalongfaultdisp_strike_plus_tplusdt_dvf_normal_minus")),
    _dalongfaultdisp_strike_minus_tplusdt_dvf_strike_plus(getMaterialPropertyByName<Real>("dalongfaultdisp_strike_minus_tplusdt_dvf_strike_plus")),
    _dalongfaultdisp_strike_minus_tplusdt_dvf_strike_minus(getMaterialPropertyByName<Real>("dalongfaultdisp_strike_minus_tplusdt_dvf_strike_minus")),
    _dalongfaultdisp_strike_minus_tplusdt_dvf_normal_plus(getMaterialPropertyByName<Real>("dalongfaultdisp_strike_minus_tplusdt_dvf_normal_plus")),
    _dalongfaultdisp_strike_minus_tplusdt_dvf_normal_minus(getMaterialPropertyByName<Real>("dalongfaultdisp_strike_minus_tplusdt_dvf_normal_minus")),
    _dalongfaultdisp_strike_plus_tplusdt_dp_plus(getMaterialPropertyByName<Real>("dalongfaultdisp_strike_plus_tplusdt_dp_plus")),
    _dalongfaultdisp_strike_plus_tplusdt_dp_minus(getMaterialPropertyByName<Real>("dalongfaultdisp_strike_plus_tplusdt_dp_minus")),
    _dalongfaultdisp_strike_minus_tplusdt_dp_plus(getMaterialPropertyByName<Real>("dalongfaultdisp_strike_minus_tplusdt_dp_plus")),
    _dalongfaultdisp_strike_minus_tplusdt_dp_minus(getMaterialPropertyByName<Real>("dalongfaultdisp_strike_minus_tplusdt_dp_minus")),
    _y_var(coupled("y_var")),
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
PoroRateStateInterfaceKernelGlobalx::computeQpResidual(Moose::DGResidualType type)
{

  Real r = 0;
  switch (type)
  {
    case Moose::Element:
      r = (_disp_strike_plus[_qp]) * _test[_i][_qp];
      break;

    case Moose::Neighbor:
      r = (_disp_strike_minus[_qp]) * _test_neighbor[_i][_qp];
      break;
  }
  return r;
  
}

Real
PoroRateStateInterfaceKernelGlobalx::computeQpJacobian(Moose::DGJacobianType type)
{
  Real jac = 0;
  switch (type)
  {
    case Moose::ElementElement:
      jac = _dalongfaultdisp_strike_plus_tplusdt_du_strike_plus[_qp] * _test[_i][_qp] * _phi[_j][_qp];
      break;
      
    case Moose::ElementNeighbor:
      jac = _dalongfaultdisp_strike_plus_tplusdt_du_strike_minus[_qp] * _test[_i][_qp] * _phi_neighbor[_j][_qp];
      break;
      
    case Moose::NeighborElement:
      jac = _dalongfaultdisp_strike_minus_tplusdt_du_strike_plus[_qp] * _test_neighbor[_i][_qp] * _phi[_j][_qp];
      break;
      
    case Moose::NeighborNeighbor:
      jac = _dalongfaultdisp_strike_minus_tplusdt_du_strike_minus[_qp] * _test_neighbor[_i][_qp] * _phi_neighbor[_j][_qp];
      break;
  }
  return jac;
}

Real
PoroRateStateInterfaceKernelGlobalx::computeQpOffDiagJacobian(Moose::DGJacobianType type, unsigned int jvar)
{
  Real jac = 0;
  
  if (jvar == _y_var)  // Coupling with y direction
  {
    switch (type)
    {
      case Moose::ElementElement:
        jac = _dalongfaultdisp_strike_plus_tplusdt_du_normal_plus[_qp] * _test[_i][_qp] * _phi[_j][_qp];
        break;
        
      case Moose::ElementNeighbor:
        jac = _dalongfaultdisp_strike_plus_tplusdt_du_normal_minus[_qp] * _test[_i][_qp] * _phi_neighbor[_j][_qp];
        break;
        
      case Moose::NeighborElement:
        jac = _dalongfaultdisp_strike_minus_tplusdt_du_normal_plus[_qp] * _test_neighbor[_i][_qp] * _phi[_j][_qp];
        break;
        
      case Moose::NeighborNeighbor:
        jac = _dalongfaultdisp_strike_minus_tplusdt_du_normal_minus[_qp] * _test_neighbor[_i][_qp] * _phi_neighbor[_j][_qp];
        break;
    }
  }
  else if (jvar == _vfx_var) // Coupling with x fluid velocity
  {
    switch (type)
    {
      case Moose::ElementElement:
        jac = _dalongfaultdisp_strike_plus_tplusdt_dvf_strike_plus[_qp] * _test[_i][_qp] * _phi_vfx[_j][_qp];
        break;
        
      case Moose::ElementNeighbor:
        jac = _dalongfaultdisp_strike_plus_tplusdt_dvf_strike_minus[_qp] * _test[_i][_qp] * _phi_vfx_neighbor[_j][_qp];
        break;
        
      case Moose::NeighborElement:
        jac = _dalongfaultdisp_strike_minus_tplusdt_dvf_strike_plus[_qp] * _test_neighbor[_i][_qp] * _phi_vfx[_j][_qp];
        break;
        
      case Moose::NeighborNeighbor:
        jac = _dalongfaultdisp_strike_minus_tplusdt_dvf_strike_minus[_qp] * _test_neighbor[_i][_qp] * _phi_vfx_neighbor[_j][_qp];
        break;
    }
  }
  else if (jvar == _vfy_var) // Coupling with y fluid velocity
  {
    switch (type)
    {
      case Moose::ElementElement:
        jac = _dalongfaultdisp_strike_plus_tplusdt_dvf_normal_plus[_qp] * _test[_i][_qp] * _phi_vfy[_j][_qp];
        break;
        
      case Moose::ElementNeighbor:
        jac = _dalongfaultdisp_strike_plus_tplusdt_dvf_normal_minus[_qp] * _test[_i][_qp] * _phi_vfy_neighbor[_j][_qp];
        break;
        
      case Moose::NeighborElement:
        jac = _dalongfaultdisp_strike_minus_tplusdt_dvf_normal_plus[_qp] * _test_neighbor[_i][_qp] * _phi_vfy[_j][_qp];
        break;
        
      case Moose::NeighborNeighbor:
        jac = _dalongfaultdisp_strike_minus_tplusdt_dvf_normal_minus[_qp] * _test_neighbor[_i][_qp] * _phi_vfy_neighbor[_j][_qp];
        break;
    }
  }
  else if (jvar == _pressure_var) // Coupling with pressure
  {
    switch (type)
    {
      case Moose::ElementElement:
        jac = _dalongfaultdisp_strike_plus_tplusdt_dp_plus[_qp] * _test[_i][_qp] * _phi_pressure[_j][_qp];
        break;
        
      case Moose::ElementNeighbor:
        jac = _dalongfaultdisp_strike_plus_tplusdt_dp_minus[_qp] * _test[_i][_qp] * _phi_pressure_neighbor[_j][_qp];
        break;
        
      case Moose::NeighborElement:
        jac = _dalongfaultdisp_strike_minus_tplusdt_dp_plus[_qp] * _test_neighbor[_i][_qp] * _phi_pressure[_j][_qp];
        break;
        
      case Moose::NeighborNeighbor:
        jac = _dalongfaultdisp_strike_minus_tplusdt_dp_minus[_qp] * _test_neighbor[_i][_qp] * _phi_pressure_neighbor[_j][_qp];
        break;
    }
  }
  return jac;
}