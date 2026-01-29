//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DGFaultSlipInterfaceKernel.h"
#include <algorithm>

registerMooseObject("farmsApp", DGFaultSlipInterfaceKernel);

InputParameters
DGFaultSlipInterfaceKernel::validParams()
{
  InputParameters params = InterfaceKernel::validParams();
  params.addClassDescription(
      "Enforces prescribed slip at fault interface using DG/Nitsche method. "
      "Makes elasticity problem linear for staggered SEAS solver.");

  params.addRequiredCoupledVar("slip_prescribed", "The prescribed slip variable");
  params.addRequiredParam<MaterialPropertyName>("shear_modulus",
                                                 "Material property for shear modulus");
  params.addParam<Real>("penalty", 1e10, "Penalty parameter for slip constraint");
  params.addParam<Real>("epsilon", 1.0, "SIPG symmetry parameter (+1 SIPG, -1 NIPG, 0 IIPG)");
  params.addParam<Real>("sigma", 6.0, "Penalty scaling factor");

  return params;
}

DGFaultSlipInterfaceKernel::DGFaultSlipInterfaceKernel(const InputParameters & parameters)
  : InterfaceKernel(parameters),
    _slip_prescribed(coupledValue("slip_prescribed")),
    _shear_modulus(getMaterialProperty<Real>("shear_modulus")),
    _penalty(getParam<Real>("penalty")),
    _epsilon(getParam<Real>("epsilon")),
    _sigma(getParam<Real>("sigma"))
{
}

Real
DGFaultSlipInterfaceKernel::computeQpResidual(Moose::DGResidualType type)
{
  // Current slip using benchmark convention: δ = u_neighbor - u_elem = u_right - u_left
  // This matches BP2 definition: δ(z,t) = u(0+,z,t) - u(0-,z,t)
  // Positive δ corresponds to right-lateral slip
  Real slip_current = _neighbor_value[_qp] - _u[_qp];

  // Slip error: δ - slip_prescribed
  Real slip_error = slip_current - _slip_prescribed[_qp];

  // Average shear modulus
  Real mu = _shear_modulus[_qp];

  // Normal vector components for gradient
  const RealVectorValue & n = _normals[_qp];

  // Average flux: {{μ ∂u/∂n}} = 0.5 * μ * (∇u_elem + ∇u_neighbor) · n
  Real avg_flux = 0.5 * mu * ((_grad_u[_qp] + _grad_neighbor_value[_qp]) * n);

  // Element size estimate (consistent with DGElasticityAntiplane)
  Real h = _current_elem_volume / _current_side_volume;

  // Penalty coefficient - use the larger of sigma*mu/h or explicit penalty
  Real kappa = std::max(_sigma * mu / h, _penalty / h);

  // SIPG formulation for enforcing [u] = slip_prescribed
  // Residual contributions similar to Dirichlet BC via Nitsche
  Real r = 0.0;

  switch (type)
  {
    case Moose::Element:
      // Element side contribution (matching MOOSE DGDiffusion sign convention):
      // -{{μ ∂u/∂n}} * v + ε * [u - slip] * (0.5*μ*∂v/∂n) + κ * [u - slip] * v
      r = -avg_flux * _test[_i][_qp];
      r += _epsilon * slip_error * 0.5 * mu * (_grad_test[_i][_qp] * n);
      r += kappa * slip_error * _test[_i][_qp];
      break;

    case Moose::Neighbor:
      // Neighbor side contribution (matching MOOSE DGDiffusion sign convention):
      // +{{μ ∂u/∂n}} * v_n + ε * [u - slip] * (0.5*μ*∂v_n/∂n) - κ * [u - slip] * v_n
      r = avg_flux * _test_neighbor[_i][_qp];
      r += _epsilon * slip_error * 0.5 * mu * (_grad_test_neighbor[_i][_qp] * n);
      r -= kappa * slip_error * _test_neighbor[_i][_qp];
      break;
  }

  return r;
}

Real
DGFaultSlipInterfaceKernel::computeQpJacobian(Moose::DGJacobianType type)
{
  Real mu = _shear_modulus[_qp];
  const RealVectorValue & n = _normals[_qp];
  Real h = _current_elem_volume / _current_side_volume;
  Real kappa = std::max(_sigma * mu / h, _penalty / h);

  Real jac = 0.0;

  // Note: slip_error = (u_neighbor - u_elem) - slip_prescribed
  // So: d(slip_error)/d(u_elem) = -1, d(slip_error)/d(u_neighbor) = +1

  switch (type)
  {
    case Moose::ElementElement:
      // d(R_elem)/d(u_elem): d(slip_error)/d(u_elem) = -1
      jac = -0.5 * mu * (_grad_phi[_j][_qp] * n) * _test[_i][_qp];
      jac += _epsilon * (-_phi[_j][_qp]) * 0.5 * mu * (_grad_test[_i][_qp] * n);
      jac += kappa * (-_phi[_j][_qp]) * _test[_i][_qp];
      break;

    case Moose::ElementNeighbor:
      // d(R_elem)/d(u_neighbor): d(slip_error)/d(u_neighbor) = +1
      jac = -0.5 * mu * (_grad_phi_neighbor[_j][_qp] * n) * _test[_i][_qp];
      jac += _epsilon * _phi_neighbor[_j][_qp] * 0.5 * mu * (_grad_test[_i][_qp] * n);
      jac += kappa * _phi_neighbor[_j][_qp] * _test[_i][_qp];
      break;

    case Moose::NeighborElement:
      // d(R_neighbor)/d(u_elem): d(slip_error)/d(u_elem) = -1
      jac = 0.5 * mu * (_grad_phi[_j][_qp] * n) * _test_neighbor[_i][_qp];
      jac += _epsilon * (-_phi[_j][_qp]) * 0.5 * mu * (_grad_test_neighbor[_i][_qp] * n);
      jac -= kappa * (-_phi[_j][_qp]) * _test_neighbor[_i][_qp];
      break;

    case Moose::NeighborNeighbor:
      // d(R_neighbor)/d(u_neighbor): d(slip_error)/d(u_neighbor) = +1
      jac = 0.5 * mu * (_grad_phi_neighbor[_j][_qp] * n) * _test_neighbor[_i][_qp];
      jac += _epsilon * _phi_neighbor[_j][_qp] * 0.5 * mu * (_grad_test_neighbor[_i][_qp] * n);
      jac -= kappa * _phi_neighbor[_j][_qp] * _test_neighbor[_i][_qp];
      break;
  }

  return jac;
}
