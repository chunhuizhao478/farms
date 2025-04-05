#include "FaultShearTractionLagrangeMultiplier.h"
registerMooseObject("farmsApp", FaultShearTractionLagrangeMultiplier);

InputParameters
FaultShearTractionLagrangeMultiplier::validParams()
{
  InputParameters params = InterfaceKernel::validParams();
  params.addRequiredCoupledVar("lambda_x", "Lagrange multiplier for shear traction");
  params.addClassDescription("Interface kernel for enforcing shear traction continuity with Lagrange multiplier");
  return params;
}

FaultShearTractionLagrangeMultiplier::FaultShearTractionLagrangeMultiplier(const InputParameters & parameters)
  : InterfaceKernel(parameters),
    _lambda_x(coupledValue("lambda_x")),
    _lambda_x_neighbor(coupledNeighborValue("lambda_x")),
    _lambda_x_var(coupled("lambda_x")),
    _lambda_x_neighbor_var(coupled("lambda_x"))
{
}

Real
FaultShearTractionLagrangeMultiplier::computeQpResidual(Moose::DGResidualType type)
{
  Real r = 0;
  switch (type)
  {
    // Residual for Element side (element) - corresponds to S𝑓+
    case Moose::DGResidualType::Element:
      // After reversing sign: +∫S𝑓+N𝑇mNpT𝑓p(t)dS
      r = _test[_i][_qp] * _lambda_x[_qp];
      break;
    // Residual for Neighbor side (neighbor) - corresponds to S𝑓-
    case Moose::DGResidualType::Neighbor:
      // After reversing sign: -∫S𝑓−N𝑇mNpT𝑓pdS
      r = -_test_neighbor[_i][_qp] * _lambda_x[_qp];
      break;
  }
  return r;
}

Real
FaultShearTractionLagrangeMultiplier::computeQpJacobian(Moose::DGJacobianType type)
{
  // The Jacobian with respect to the displacement variables is zero
  // since the Lagrange multiplier terms don't directly depend on displacement
  Real jac = 0;
  switch (type)
  {
    case Moose::DGJacobianType::ElementElement:
      jac = _test[_i][_qp] * _phi[_j][_qp]* _lambda_x[_qp];
      break;
    case Moose::DGJacobianType::NeighborNeighbor:
      jac = - _test_neighbor[_i][_qp] * _phi_neighbor[_j][_qp]* _lambda_x[_qp];
      break;
    case Moose::DGJacobianType::ElementNeighbor:
      jac = _test[_i][_qp] * _phi_neighbor[_j][_qp]* _lambda_x[_qp];
      break;
    case Moose::DGJacobianType::NeighborElement:
      jac = - _test_neighbor[_i][_qp] * _phi[_j][_qp]* _lambda_x[_qp];
      break;
  }
  return jac;
}

Real
FaultShearTractionLagrangeMultiplier::computeQpOffDiagJacobian(Moose::DGJacobianType type, unsigned int jvar)
{
  // Jacobian contribution from Lagrange multiplier
  if (jvar == _lambda_x_var)
  {
    switch (type)
    {
      case Moose::DGJacobianType::ElementElement:
        return _test[_i][_qp] * _phi[_j][_qp];
      case Moose::DGJacobianType::NeighborNeighbor:
        return -_test_neighbor[_i][_qp] * _phi_neighbor[_j][_qp]; // No contribution from element lambda to neighbor
      default:
        return 0.0;
    }
  }
   return 0.0;
}