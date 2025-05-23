#include "FaultSlipInterfaceKernel.h"
registerMooseObject("farmsApp", FaultSlipInterfaceKernel);

InputParameters
FaultSlipInterfaceKernel::validParams()
{
  InputParameters params = InterfaceKernel::validParams();
  params.addClassDescription("Interface kernel to enforce fault slip condition with penalty");
  params.addRequiredCoupledVar("slip", "The slip value to enforce at the interface");
  params.addRequiredCoupledVar("coupled_disp", "Variable on the element side");

  return params;
}

FaultSlipInterfaceKernel::FaultSlipInterfaceKernel(const InputParameters & parameters)
  : InterfaceKernel(parameters),
    _slip(coupledValue("slip")),
    _slip_var(coupled("slip")),
    _coupled_disp(coupledValue("coupled_disp")),
    _coupled_disp_num(coupled("coupled_disp")),
    _coupled_disp_neighbor(coupledNeighborValue("coupled_disp")),
    _coupled_disp_neighbor_num(coupled("coupled_disp"))
{
}

Real
FaultSlipInterfaceKernel::computeQpResidual(Moose::DGResidualType type)
{
  // Calculate the jump in displacement
  Real jump = _coupled_disp[_qp] - _coupled_disp_neighbor[_qp];
  
  // Residual = (jump - target_slip)
  Real r =  (jump - _slip[_qp]);
  
  switch (type)
  {
    case Moose::DGResidualType::Element:
      return -_test[_i][_qp] * r;
    
    case Moose::DGResidualType::Neighbor:
      return _test_neighbor[_i][_qp] * r;
    
    default:
      return 0.0;
  }
}

Real
FaultSlipInterfaceKernel::computeQpJacobian(Moose::DGJacobianType type)
{
  // Only contribute to Jacobian if this variable is the displacement
  if (_var.number() != _coupled_disp_num)
    return 0.0;
  // Calculate the jump in displacement
  Real jump = _coupled_disp[_qp] - _coupled_disp_neighbor[_qp];  
  switch (type)
  {
    case Moose::DGJacobianType::ElementElement:
      return _test[_i][_qp] *  _phi[_j][_qp] * jump;
    
    case Moose::DGJacobianType::NeighborNeighbor:
      return -_test_neighbor[_i][_qp] *  _phi_neighbor[_j][_qp]* jump;
    
    case Moose::DGJacobianType::ElementNeighbor:
      return -_test[_i][_qp] *  _phi_neighbor[_j][_qp]* jump;
    
    case Moose::DGJacobianType::NeighborElement:
      return _test_neighbor[_i][_qp] *  _phi[_j][_qp]* jump;
    
    default:
      return 0.0;
  }
}

Real
FaultSlipInterfaceKernel::computeQpOffDiagJacobian(Moose::DGJacobianType type, unsigned int jvar)
{
  
  if (jvar == _coupled_disp_num)
  {
    // Jacobian for coupled displacement
    switch (type)
    {
      case Moose::DGJacobianType::ElementElement:
        return _test[_i][_qp] *  _phi[_j][_qp];
      
      case Moose::DGJacobianType::NeighborNeighbor:
        return -_test_neighbor[_i][_qp] *  _phi_neighbor[_j][_qp];
      
      case Moose::DGJacobianType::ElementNeighbor:
        return -_test[_i][_qp] *  _phi_neighbor[_j][_qp];
      
      case Moose::DGJacobianType::NeighborElement:
        return _test_neighbor[_i][_qp] *  _phi[_j][_qp];
      
      default:
        return 0.0;
    }
  }
  
  return 0.0;
}