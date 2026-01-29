//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DGAntiplaneFaultInterfaceKernel.h"

registerMooseObject("farmsApp", DGAntiplaneFaultInterfaceKernel);

InputParameters
DGAntiplaneFaultInterfaceKernel::validParams()
{
  InputParameters params = InterfaceKernel::validParams();
  params.addClassDescription(
      "Fault interface kernel for DG antiplane shear elasticity in SEAS simulations. "
      "Applies fault traction from rate-state friction law and allows discontinuous slip.");

  params.addRequiredParam<MaterialPropertyName>(
      "fault_traction", "Material property for fault shear traction from friction law");

  params.addParam<MaterialPropertyName>(
      "dtraction_dslip", "Material property for derivative of traction w.r.t. slip");
  params.addParam<MaterialPropertyName>(
      "dtraction_dslip_rate", "Material property for derivative of traction w.r.t. slip rate");

  params.addParam<bool>("use_penalty", false,
                        "Whether to use penalty enforcement for slip-traction coupling");
  params.addParam<Real>("penalty", 1e10, "Penalty parameter for slip enforcement");

  params.addParam<MaterialPropertyName>("shear_modulus", "shear_modulus",
                                        "Shear modulus for penalty calculation");

  return params;
}

DGAntiplaneFaultInterfaceKernel::DGAntiplaneFaultInterfaceKernel(const InputParameters & parameters)
  : InterfaceKernel(parameters),
    _fault_traction(getMaterialProperty<Real>("fault_traction")),
    _dtraction_dslip(isParamValid("dtraction_dslip")
                         ? &getMaterialProperty<Real>("dtraction_dslip")
                         : nullptr),
    _dtraction_dslip_rate(isParamValid("dtraction_dslip_rate")
                              ? &getMaterialProperty<Real>("dtraction_dslip_rate")
                              : nullptr),
    _use_penalty(getParam<bool>("use_penalty")),
    _penalty(getParam<Real>("penalty")),
    _mu(_use_penalty ? &getMaterialProperty<Real>("shear_modulus") : nullptr),
    _mu_neighbor(_use_penalty ? &getNeighborMaterialProperty<Real>("shear_modulus") : nullptr)
{
}

Real
DGAntiplaneFaultInterfaceKernel::computeQpResidual(Moose::DGResidualType type)
{
  Real r = 0.0;

  // Get fault traction from friction law (computed by material)
  const Real tau = _fault_traction[_qp];

  // The weak form for fault interface is:
  // ∫_Γf τ * [[v]] dA = ∫_Γf τ * (v+ - v-) dA
  //
  // This means:
  // - On Element side (+): residual = τ * v+
  // - On Neighbor side (-): residual = -τ * v-
  //
  // The traction τ acts in the direction that opposes relative motion.
  // For antiplane shear, positive τ means traction in +z direction.

  switch (type)
  {
    case Moose::Element:
      // Residual on element side: τ * v_elem
      // This comes from: -∫ σ·n * v dA = -∫ τ * v dA (moving to RHS gives +τ*v)
      // Note: Sign convention - traction acts on the + side
      r = tau * _test[_i][_qp];
      break;

    case Moose::Neighbor:
      // Residual on neighbor side: -τ * v_neighbor
      // By Newton's 3rd law, traction on - side is opposite
      r = -tau * _test_neighbor[_i][_qp];
      break;
  }

  // Optional: Add penalty term to enforce slip-traction coupling
  // This helps with convergence when slip and traction are strongly coupled
  if (_use_penalty)
  {
    // Jump in displacement (slip)
    const Real slip = _u[_qp] - _neighbor_value[_qp];

    // Penalty contribution: penalty * (u+ - u-) * [[v]]
    switch (type)
    {
      case Moose::Element:
        r += _penalty * slip * _test[_i][_qp];
        break;
      case Moose::Neighbor:
        r -= _penalty * slip * _test_neighbor[_i][_qp];
        break;
    }
  }

  return r;
}

Real
DGAntiplaneFaultInterfaceKernel::computeQpJacobian(Moose::DGJacobianType type)
{
  Real r = 0.0;

  // Jacobian contributions from traction dependence on displacement
  // τ = τ(s, V) where s = [[u]] is slip and V = ds/dt is slip rate
  //
  // dτ/du = dτ/ds * ds/du = dτ/ds (since s = u+ - u-)
  // On Element side: ds/du+ = 1
  // On Neighbor side: ds/du- = -1

  if (_dtraction_dslip)
  {
    const Real dtau_dslip = (*_dtraction_dslip)[_qp];

    switch (type)
    {
      case Moose::ElementElement:
        // d(τ * v+)/du+ = dτ/ds * ds/du+ * v+ = dτ/ds * 1 * v+ * phi_j
        r = dtau_dslip * _phi[_j][_qp] * _test[_i][_qp];
        break;

      case Moose::ElementNeighbor:
        // d(τ * v+)/du- = dτ/ds * ds/du- * v+ = dτ/ds * (-1) * v+ * phi_neighbor_j
        r = -dtau_dslip * _phi_neighbor[_j][_qp] * _test[_i][_qp];
        break;

      case Moose::NeighborElement:
        // d(-τ * v-)/du+ = -dτ/ds * ds/du+ * v- = -dτ/ds * 1 * v- * phi_j
        r = -dtau_dslip * _phi[_j][_qp] * _test_neighbor[_i][_qp];
        break;

      case Moose::NeighborNeighbor:
        // d(-τ * v-)/du- = -dτ/ds * ds/du- * v- = -dτ/ds * (-1) * v- * phi_neighbor_j
        r = dtau_dslip * _phi_neighbor[_j][_qp] * _test_neighbor[_i][_qp];
        break;
    }
  }

  // Penalty term Jacobian
  if (_use_penalty)
  {
    switch (type)
    {
      case Moose::ElementElement:
        r += _penalty * _phi[_j][_qp] * _test[_i][_qp];
        break;

      case Moose::ElementNeighbor:
        r -= _penalty * _phi_neighbor[_j][_qp] * _test[_i][_qp];
        break;

      case Moose::NeighborElement:
        r -= _penalty * _phi[_j][_qp] * _test_neighbor[_i][_qp];
        break;

      case Moose::NeighborNeighbor:
        r += _penalty * _phi_neighbor[_j][_qp] * _test_neighbor[_i][_qp];
        break;
    }
  }

  return r;
}
