//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DGElasticTractionInterfaceUO.h"
#include <algorithm>

registerMooseObject("farmsApp", DGElasticTractionInterfaceUO);

InputParameters
DGElasticTractionInterfaceUO::validParams()
{
  InputParameters params = InterfaceQpUserObjectBase::validParams();
  params.addClassDescription(
      "Computes DG-consistent elastic traction at fault interfaces for SEAS problems. "
      "Uses averaged gradient and penalty term following Tandem's SIPG formulation.");

  params.addRequiredCoupledVar("displacement", "The displacement variable");
  params.addCoupledVar("slip_prescribed", "The prescribed slip AuxVariable (optional)");
  params.addRequiredParam<Real>("shear_modulus", "Shear modulus mu (Pa)");
  params.addParam<Real>("sigma", 6.0, "DG penalty scaling factor");
  params.addParam<Real>("penalty", 1e10, "Explicit penalty parameter");
  params.addParam<Real>("tau_pre", 0.0, "Pre-stress/initial traction to add (Pa)");

  return params;
}

DGElasticTractionInterfaceUO::DGElasticTractionInterfaceUO(const InputParameters & parameters)
  : InterfaceQpUserObjectBase(parameters),
    _u(coupledValue("displacement")),
    _u_neighbor(coupledNeighborValue("displacement")),
    _grad_u(coupledGradient("displacement")),
    _grad_u_neighbor(coupledNeighborGradient("displacement")),
    _slip_prescribed(isParamValid("slip_prescribed") ? &coupledValue("slip_prescribed") : nullptr),
    _shear_modulus(getParam<Real>("shear_modulus")),
    _sigma(getParam<Real>("sigma")),
    _penalty(getParam<Real>("penalty")),
    _tau_pre(getParam<Real>("tau_pre")),
    _assembly(_subproblem.assembly(_tid, 0))
{
  if (_shear_modulus <= 0.0)
    mooseError("Shear modulus must be positive");
}

Real
DGElasticTractionInterfaceUO::computeRealValue(const unsigned int qp)
{
  // DG elastic traction following Tandem's Poisson operator:
  // traction = mu * {{grad_u}} . n - penalty * slip_error
  //
  // Physical reasoning:
  // - If slip_error > 0 (more slip than expected): stress is released, traction decreases
  // - If slip_error < 0 (less slip, fault locked): stress builds up, traction increases

  // Get the interface normal (points from elem to neighbor)
  const RealVectorValue & n = _normals[qp];

  // Average gradient projected onto normal: tau = mu * {{grad_u}} . n
  Real avg_grad_dot_n = 0.5 * ((_grad_u[qp] + _grad_u_neighbor[qp]) * n);
  Real avg_stress_traction = _shear_modulus * avg_grad_dot_n;

  // Compute slip error for penalty correction (if slip is prescribed)
  Real slip_error = 0.0;
  if (_slip_prescribed)
  {
    // Using MOOSE convention: slip = u_neighbor - u_elem
    Real slip_current = _u_neighbor[qp] - _u[qp];
    slip_error = slip_current - (*_slip_prescribed)[qp];
  }

  // Element size estimate for penalty computation
  Real h = _current_elem_volume / _current_side_volume;

  // Penalty coefficient (same formula as DGFaultSlipInterfaceKernel)
  Real penalty_coef = std::max(_sigma * _shear_modulus / h, _penalty / h);

  // DG-consistent traction formula following Tandem
  // traction = avg_stress_traction - penalty * slip_error + tau_pre
  Real traction = avg_stress_traction - penalty_coef * slip_error + _tau_pre;

  return traction;
}
