#include "FaultNormalTractionLagrangeMultiplier.h"
registerMooseObject("farmsApp", FaultNormalTractionLagrangeMultiplier);

InputParameters
FaultNormalTractionLagrangeMultiplier::validParams()
{
  InputParameters params = InterfaceKernel::validParams();
  params.addRequiredCoupledVar("lambda_y", "Lagrange multiplier for normal direction");
  params.addClassDescription("Interface kernel for enforcing normal traction continuity with Lagrange multiplier");
  return params;
}

FaultNormalTractionLagrangeMultiplier::FaultNormalTractionLagrangeMultiplier(const InputParameters & parameters)
  : InterfaceKernel(parameters),
    _lambda_y(coupledValue("lambda_y")),
    _lambda_y_neighbor(coupledNeighborValue("lambda_y")),
    _lambda_y_var(coupled("lambda_y")),
    _lambda_y_neighbor_var(coupled("lambda_y"))
{
}

Real
FaultNormalTractionLagrangeMultiplier::computeQpResidual(Moose::DGResidualType type)
{
  Real r = 0;
  switch (type)
  {
    // Residual for Element side (element) - corresponds to S𝑓+
    case Moose::DGResidualType::Element:
      // After reversing sign: +∫S𝑓+N𝑇mNpT𝑓p(t)dS
      r = _test[_i][_qp] * _lambda_y[_qp];
      break;
    // Residual for Neighbor side (neighbor) - corresponds to S𝑓-
    case Moose::DGResidualType::Neighbor:
      // After reversing sign: -∫S𝑓−N𝑇mNpT𝑓pdS
      // Note: For normal traction, we use a positive sign for neighbor residual
      r = _test_neighbor[_i][_qp] * _lambda_y[_qp];
      break;
  }
  return r;
}

Real
FaultNormalTractionLagrangeMultiplier::computeQpJacobian(Moose::DGJacobianType type)
{
  // The Jacobian with respect to the displacement variables is zero
  // since the Lagrange multiplier terms don't directly depend on displacement
   Real jac = 0;
  switch (type)
  {
    case Moose::DGJacobianType::ElementElement:
      jac = _test[_i][_qp] * _phi[_j][_qp]* _lambda_y[_qp];
      break;
    case Moose::DGJacobianType::NeighborNeighbor:
      jac = - _test_neighbor[_i][_qp] * _phi_neighbor[_j][_qp]* _lambda_y[_qp];
      break;
    case Moose::DGJacobianType::ElementNeighbor:
      jac = _test[_i][_qp] * _phi_neighbor[_j][_qp]* _lambda_y[_qp];
      break;
    case Moose::DGJacobianType::NeighborElement:
      jac = - _test_neighbor[_i][_qp] * _phi[_j][_qp]* _lambda_y[_qp];
      break;
  }
  return jac;
}

Real
FaultNormalTractionLagrangeMultiplier::computeQpOffDiagJacobian(Moose::DGJacobianType type, unsigned int jvar)
{
  // Jacobian contribution from Lagrange multiplier
  if (jvar == _lambda_y_var)
  {
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
      
    }
  }
  return 0.0;
}