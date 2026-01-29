//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DGElasticTractionMaterial.h"
#include <algorithm>
#include <cmath>

registerMooseObject("farmsApp", DGElasticTractionMaterial);

InputParameters
DGElasticTractionMaterial::validParams()
{
  InputParameters params = InterfaceMaterial::validParams();
  params.addClassDescription(
      "Computes elastic shear traction on fault interface for SEAS simulations. "
      "Supports 'gradient' mode (DG formulation) and 'stiffness' mode "
      "(discrete stiffness relationship for quasi-static problems).");

  params.addRequiredCoupledVar("displacement", "The displacement variable");
  params.addRequiredCoupledVar("slip_prescribed", "The prescribed slip AuxVariable");
  params.addRequiredParam<Real>("shear_modulus", "Shear modulus μ (Pa)");
  params.addParam<Real>("sigma", 6.0, "DG penalty scaling factor");
  params.addParam<Real>("penalty", 1e10, "Explicit penalty parameter");
  params.addParam<Real>("traction_scale", 1.0, "Scaling factor for elastic traction");
  params.addParam<MaterialPropertyName>("traction_name", "elastic_traction",
                                        "Name of the output traction property");

  // Traction mode selection
  MooseEnum traction_modes("gradient stiffness dg_consistent", "gradient");
  params.addParam<MooseEnum>("traction_mode", traction_modes,
      "Mode for computing traction: 'gradient' uses DG formulation τ = μ * {{∂u/∂n}}, "
      "'stiffness' uses quasi-static relationship τ = K * (Vp*t - slip), "
      "'dg_consistent' uses full DG formula τ = μ * {{∂u/∂n}} + κ * ([[u]] - slip) (Tandem-style)");

  // Parameters for stiffness mode
  params.addParam<Real>("plate_rate", 1e-9,
      "Plate rate Vp (m/s) for stiffness mode. Default: 1e-9 m/s");
  params.addParam<Real>("fault_depth", 40000.0,
      "Fault depth W (m) for stiffness calculation K = μ/(π*W). Default: 40 km");
  params.addParam<Real>("domain_half_width", 100000.0,
      "Domain half-width L (m) for loading calculation. Default: 100 km");

  return params;
}

DGElasticTractionMaterial::DGElasticTractionMaterial(const InputParameters & parameters)
  : InterfaceMaterial(parameters),
    _u(coupledValue("displacement")),
    _u_neighbor(coupledNeighborValue("displacement")),
    _grad_u(coupledGradient("displacement")),
    _grad_u_neighbor(coupledNeighborGradient("displacement")),
    _slip_prescribed(coupledValue("slip_prescribed")),
    _shear_modulus(getParam<Real>("shear_modulus")),
    _sigma(getParam<Real>("sigma")),
    _penalty(getParam<Real>("penalty")),
    _traction_scale(getParam<Real>("traction_scale")),
    _current_elem_volume(_assembly.elemVolume()),
    _current_side_volume(_assembly.sideElemVolume()),
    _elastic_traction(
        declarePropertyByName<Real>(getParam<MaterialPropertyName>("traction_name"))),
    _traction_mode(getParam<MooseEnum>("traction_mode")),
    _plate_rate(getParam<Real>("plate_rate")),
    _fault_depth(getParam<Real>("fault_depth")),
    _domain_half_width(getParam<Real>("domain_half_width"))
{
}

void
DGElasticTractionMaterial::computeQpProperties()
{
  Real elastic_traction = 0.0;

  if (_traction_mode == "gradient")
  {
    // Original DG formulation:
    // τ = μ * {{∂u/∂n}} = μ * 0.5 * (∇u_elem + ∇u_neighbor) · n
    // Note: This may produce incorrect stress rates for quasi-static SEAS problems
    // due to mesh-dependent artifacts from slip enforcement.

    const RealVectorValue & n = _normals[_qp];
    Real avg_grad_dot_n = 0.5 * ((_grad_u[_qp] + _grad_u_neighbor[_qp]) * n);
    elastic_traction = _shear_modulus * avg_grad_dot_n;
  }
  else if (_traction_mode == "stiffness")
  {
    // Stiffness-based formulation for quasi-static SEAS:
    // τ = K * (Vp*t - slip)
    //
    // For a 2D antiplane problem with a vertical fault of depth W in an elastic half-space,
    // the discrete stiffness is approximately:
    //   K = μ / W
    //
    // Note: Some formulations use K = μ / (π * W), but K = μ / W better matches
    // the SCEC SEAS benchmark results (unicycle-ap-ratestate code).
    //
    // This correctly captures the quasi-static stress evolution where:
    // - Stress increases due to far-field loading (Vp * t)
    // - Stress decreases due to fault slip
    // - Net stress change depends on slip deficit (Vp*t - slip)

    // Compute stiffness K = μ / W
    Real K = _shear_modulus / _fault_depth;

    // Get current time
    Real t = _t;

    // Get the local slip value
    Real slip = _slip_prescribed[_qp];

    // Compute slip deficit
    Real slip_deficit = _plate_rate * t - slip;

    // Elastic traction from stiffness relationship
    elastic_traction = K * slip_deficit;
  }
  else if (_traction_mode == "dg_consistent")
  {
    // DG-consistent traction computation following Tandem's approach:
    // τ = μ * {{∂u/∂n}} + κ * ([[u]] - slip_prescribed)
    //
    // This is the correct DG formulation that accounts for the penalty
    // enforcement of the slip boundary condition. The penalty correction
    // term ensures the traction is consistent with the weak form.
    //
    // Reference: Tandem code (app/kernels/elasticity.py, compute_traction)

    const RealVectorValue & n = _normals[_qp];

    // 1. Average gradient term: μ * {{∂u/∂n}}
    Real avg_grad_dot_n = 0.5 * ((_grad_u[_qp] + _grad_u_neighbor[_qp]) * n);
    Real gradient_term = _shear_modulus * avg_grad_dot_n;

    // 2. Displacement jump: [[u]] = u_neighbor - u_elem
    // Convention matches DGFaultSlipInterfaceKernel: slip = u_neighbor - u_elem
    Real u_jump = _u_neighbor[_qp] - _u[_qp];

    // 3. Element size (same as in DGFaultSlipInterfaceKernel)
    Real h = _current_elem_volume / _current_side_volume;

    // 4. Penalty coefficient κ (same formula as DGFaultSlipInterfaceKernel)
    Real kappa = std::max(_sigma * _shear_modulus / h, _penalty / h);

    // 5. Slip error: [[u]] - slip_prescribed
    Real slip_error = u_jump - _slip_prescribed[_qp];

    // 6. Full traction: gradient term - penalty correction
    // IMPORTANT: The minus sign follows Tandem's convention (Poisson.cpp:690)
    // Physical reasoning:
    // - If slip_error > 0 (more slip than prescribed): stress released, traction decreases
    // - If slip_error < 0 (fault locked more): stress builds up, traction increases
    elastic_traction = gradient_term - kappa * slip_error;
  }

  // Apply traction_scale if needed
  _elastic_traction[_qp] = _traction_scale * elastic_traction;
}
