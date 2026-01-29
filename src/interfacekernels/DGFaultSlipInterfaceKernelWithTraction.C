//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DGFaultSlipInterfaceKernelWithTraction.h"
#include "MooseVariable.h"
#include <algorithm>

registerMooseObject("farmsApp", DGFaultSlipInterfaceKernelWithTraction);

InputParameters
DGFaultSlipInterfaceKernelWithTraction::validParams()
{
  InputParameters params = InterfaceKernel::validParams();
  params.addClassDescription(
      "Enforces prescribed slip at fault interface using DG/Nitsche method "
      "and outputs the computed traction to an AuxVariable.");

  params.addRequiredCoupledVar("slip_prescribed", "The prescribed slip variable");
  params.addRequiredParam<MaterialPropertyName>("shear_modulus",
                                                 "Material property for shear modulus");
  params.addParam<Real>("penalty", 1e10, "Penalty parameter for slip constraint");
  params.addParam<Real>("epsilon", 1.0, "SIPG symmetry parameter (+1 SIPG, -1 NIPG, 0 IIPG)");
  params.addParam<Real>("sigma", 6.0, "Penalty scaling factor");
  params.addParam<Real>("tau_pre", 0.0, "Pre-stress to add to computed traction (Pa)");
  params.addRequiredCoupledVar("traction_var", "AuxVariable to write the computed traction to");

  return params;
}

DGFaultSlipInterfaceKernelWithTraction::DGFaultSlipInterfaceKernelWithTraction(
    const InputParameters & parameters)
  : InterfaceKernel(parameters),
    _slip_prescribed(coupledValue("slip_prescribed")),
    _shear_modulus(getMaterialProperty<Real>("shear_modulus")),
    _penalty(getParam<Real>("penalty")),
    _epsilon(getParam<Real>("epsilon")),
    _sigma(getParam<Real>("sigma")),
    _tau_pre(getParam<Real>("tau_pre")),
    _traction_var(getVar("traction_var", 0))
{
  if (!_traction_var)
    mooseError("traction_var must be specified");
}

Real
DGFaultSlipInterfaceKernelWithTraction::computeTraction() const
{
  // Compute traction using the same values as the residual computation
  //
  // Following Tandem's convention (Poisson.cpp traction_skeleton):
  //   τ = μ * {{∂u/∂n}} - κ * ([[u]] - slip)
  //
  // Physical interpretation:
  // - {{∂u/∂n}} captures the average stress from the displacement gradient
  // - The penalty term corrects for the slip constraint enforcement
  // - Minus sign: if [[u]] > slip (more slip than prescribed), traction decreases

  Real mu = _shear_modulus[_qp];
  const RealVectorValue & n = _normals[_qp];

  // Average gradient dotted with normal: {{∂u/∂n}}
  Real avg_grad_dot_n = 0.5 * ((_grad_u[_qp] + _grad_neighbor_value[_qp]) * n);

  // Gradient term: μ * {{∂u/∂n}}
  Real gradient_term = mu * avg_grad_dot_n;

  // Slip and slip error
  Real slip_current = _neighbor_value[_qp] - _u[_qp];  // [[u]] = u_neighbor - u_elem
  Real slip_error = slip_current - _slip_prescribed[_qp];

  // Element size and penalty coefficient (same as residual computation)
  Real h = _current_elem_volume / _current_side_volume;
  Real kappa = std::max(_sigma * mu / h, _penalty / h);

  // Traction with penalty correction (MINUS sign following Tandem)
  Real traction = gradient_term - kappa * slip_error + _tau_pre;

  return traction;
}

Real
DGFaultSlipInterfaceKernelWithTraction::computeQpResidual(Moose::DGResidualType type)
{
  // Compute and store traction (only once per qp, on element side)
  if (type == Moose::Element && _i == 0)
  {
    Real traction = computeTraction();
    // Write to the AuxVariable at this quadrature point
    // Using setNodalValue for CONSTANT MONOMIAL variables on the boundary
    _traction_var->setNodalValue(traction);
  }

  // Standard slip BC residual computation (same as DGFaultSlipInterfaceKernel)
  Real slip_current = _neighbor_value[_qp] - _u[_qp];
  Real slip_error = slip_current - _slip_prescribed[_qp];
  Real mu = _shear_modulus[_qp];
  const RealVectorValue & n = _normals[_qp];
  Real avg_flux = 0.5 * mu * ((_grad_u[_qp] + _grad_neighbor_value[_qp]) * n);
  Real h = _current_elem_volume / _current_side_volume;
  Real kappa = std::max(_sigma * mu / h, _penalty / h);

  Real r = 0.0;

  switch (type)
  {
    case Moose::Element:
      r = -avg_flux * _test[_i][_qp];
      r += _epsilon * slip_error * 0.5 * mu * (_grad_test[_i][_qp] * n);
      r += kappa * slip_error * _test[_i][_qp];
      break;

    case Moose::Neighbor:
      r = avg_flux * _test_neighbor[_i][_qp];
      r += _epsilon * slip_error * 0.5 * mu * (_grad_test_neighbor[_i][_qp] * n);
      r -= kappa * slip_error * _test_neighbor[_i][_qp];
      break;
  }

  return r;
}

Real
DGFaultSlipInterfaceKernelWithTraction::computeQpJacobian(Moose::DGJacobianType type)
{
  Real mu = _shear_modulus[_qp];
  const RealVectorValue & n = _normals[_qp];
  Real h = _current_elem_volume / _current_side_volume;
  Real kappa = std::max(_sigma * mu / h, _penalty / h);

  Real jac = 0.0;

  switch (type)
  {
    case Moose::ElementElement:
      jac = -0.5 * mu * (_grad_phi[_j][_qp] * n) * _test[_i][_qp];
      jac += _epsilon * (-_phi[_j][_qp]) * 0.5 * mu * (_grad_test[_i][_qp] * n);
      jac += kappa * (-_phi[_j][_qp]) * _test[_i][_qp];
      break;

    case Moose::ElementNeighbor:
      jac = -0.5 * mu * (_grad_phi_neighbor[_j][_qp] * n) * _test[_i][_qp];
      jac += _epsilon * _phi_neighbor[_j][_qp] * 0.5 * mu * (_grad_test[_i][_qp] * n);
      jac += kappa * _phi_neighbor[_j][_qp] * _test[_i][_qp];
      break;

    case Moose::NeighborElement:
      jac = 0.5 * mu * (_grad_phi[_j][_qp] * n) * _test_neighbor[_i][_qp];
      jac += _epsilon * (-_phi[_j][_qp]) * 0.5 * mu * (_grad_test_neighbor[_i][_qp] * n);
      jac -= kappa * (-_phi[_j][_qp]) * _test_neighbor[_i][_qp];
      break;

    case Moose::NeighborNeighbor:
      jac = 0.5 * mu * (_grad_phi_neighbor[_j][_qp] * n) * _test_neighbor[_i][_qp];
      jac += _epsilon * _phi_neighbor[_j][_qp] * 0.5 * mu * (_grad_test_neighbor[_i][_qp] * n);
      jac -= kappa * _phi_neighbor[_j][_qp] * _test_neighbor[_i][_qp];
      break;
  }

  return jac;
}
