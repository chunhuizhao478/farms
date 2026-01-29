//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DGElasticityAntiplane.h"
#include "MooseVariableFE.h"
#include "libmesh/utility.h"

registerMooseObject("farmsApp", DGElasticityAntiplane);

InputParameters
DGElasticityAntiplane::validParams()
{
  InputParameters params = DGKernel::validParams();
  params.addClassDescription(
      "Implements SIPG (Symmetric Interior Penalty Galerkin) method for 2D antiplane "
      "shear elasticity. Based on the Tandem paper formulation for SEAS simulations.");
  params.addParam<Real>("epsilon", 1.0,
                        "Symmetry parameter: 1 = SIPG (symmetric), -1 = NIPG (non-symmetric), "
                        "0 = IIPG (incomplete)");
  params.addParam<Real>("sigma", 6.0, "Penalty parameter multiplier for stability");
  params.addParam<MaterialPropertyName>("shear_modulus", "shear_modulus",
                                        "The shear modulus material property name");
  return params;
}

DGElasticityAntiplane::DGElasticityAntiplane(const InputParameters & parameters)
  : DGKernel(parameters),
    _epsilon(getParam<Real>("epsilon")),
    _sigma(getParam<Real>("sigma")),
    _mu(getMaterialProperty<Real>("shear_modulus")),
    _mu_neighbor(getNeighborMaterialProperty<Real>("shear_modulus"))
{
}

Real
DGElasticityAntiplane::computeQpResidual(Moose::DGResidualType type)
{
  Real r = 0.0;

  // Compute characteristic element size h
  // h = element_volume / face_area * 1/p^2 (p is polynomial order)
  const int elem_b_order = std::max(libMesh::Order(1), _var.order());
  const Real h_elem =
      _current_elem_volume / _current_side_volume * 1.0 / Utility::pow<2>(elem_b_order);

  // Average shear modulus
  const Real mu_avg = 0.5 * (_mu[_qp] + _mu_neighbor[_qp]);

  // Jump in solution: [[u]] = u_elem - u_neighbor
  const Real jump_u = _u[_qp] - _u_neighbor[_qp];

  // Average of flux: {{mu * grad(u) . n}}
  const Real avg_flux = 0.5 * (_mu[_qp] * _grad_u[_qp] * _normals[_qp] +
                               _mu_neighbor[_qp] * _grad_u_neighbor[_qp] * _normals[_qp]);

  switch (type)
  {
    case Moose::Element:
      // -{{mu*grad(u).n}} * v_elem  (consistency term)
      r -= avg_flux * _test[_i][_qp];

      // +epsilon * [[u]] * {{mu*grad(v).n}}  (symmetry term, matches MOOSE DGDiffusion)
      // For element side, {{mu*grad(v).n}} = 0.5 * mu_elem * grad(v_elem).n
      r += _epsilon * 0.5 * jump_u * _mu[_qp] * (_grad_test[_i][_qp] * _normals[_qp]);

      // sigma/h * [[u]] * v_elem  (penalty term)
      r += _sigma * mu_avg / h_elem * jump_u * _test[_i][_qp];
      break;

    case Moose::Neighbor:
      // +{{mu*grad(u).n}} * v_neighbor  (consistency term, opposite sign)
      r += avg_flux * _test_neighbor[_i][_qp];

      // +epsilon * [[u]] * {{mu*grad(v).n}}  (symmetry term, matches MOOSE DGDiffusion)
      // For neighbor side, {{mu*grad(v).n}} = 0.5 * mu_neighbor * grad(v_neighbor).n
      r += _epsilon * 0.5 * jump_u * _mu_neighbor[_qp] *
           (_grad_test_neighbor[_i][_qp] * _normals[_qp]);

      // -sigma/h * [[u]] * v_neighbor  (penalty term, opposite sign)
      r -= _sigma * mu_avg / h_elem * jump_u * _test_neighbor[_i][_qp];
      break;
  }

  return r;
}

Real
DGElasticityAntiplane::computeQpJacobian(Moose::DGJacobianType type)
{
  Real r = 0.0;

  // Compute characteristic element size h
  const int elem_b_order = std::max(libMesh::Order(1), _var.order());
  const Real h_elem =
      _current_elem_volume / _current_side_volume * 1.0 / Utility::pow<2>(elem_b_order);

  // Average shear modulus
  const Real mu_avg = 0.5 * (_mu[_qp] + _mu_neighbor[_qp]);

  switch (type)
  {
    case Moose::ElementElement:
      // d/du_j of: -0.5 * mu_elem * grad(u).n * v_i
      r -= 0.5 * _mu[_qp] * (_grad_phi[_j][_qp] * _normals[_qp]) * _test[_i][_qp];

      // d/du_j of: +epsilon * 0.5 * (u - u_neighbor) * mu_elem * grad(v_i).n
      r += _epsilon * 0.5 * _phi[_j][_qp] * _mu[_qp] * (_grad_test[_i][_qp] * _normals[_qp]);

      // d/du_j of: sigma/h * (u - u_neighbor) * v_i
      r += _sigma * mu_avg / h_elem * _phi[_j][_qp] * _test[_i][_qp];
      break;

    case Moose::ElementNeighbor:
      // d/du_neighbor_j of: -0.5 * mu_neighbor * grad(u_neighbor).n * v_i
      r -= 0.5 * _mu_neighbor[_qp] * (_grad_phi_neighbor[_j][_qp] * _normals[_qp]) * _test[_i][_qp];

      // d/du_neighbor_j of: +epsilon * 0.5 * (u - u_neighbor) * mu_elem * grad(v_i).n
      r += _epsilon * 0.5 * (-_phi_neighbor[_j][_qp]) * _mu[_qp] *
           (_grad_test[_i][_qp] * _normals[_qp]);

      // d/du_neighbor_j of: sigma/h * (u - u_neighbor) * v_i
      r += _sigma * mu_avg / h_elem * (-_phi_neighbor[_j][_qp]) * _test[_i][_qp];
      break;

    case Moose::NeighborElement:
      // d/du_j of: +0.5 * mu_elem * grad(u).n * v_neighbor_i
      r += 0.5 * _mu[_qp] * (_grad_phi[_j][_qp] * _normals[_qp]) * _test_neighbor[_i][_qp];

      // d/du_j of: +epsilon * 0.5 * (u - u_neighbor) * mu_neighbor * grad(v_neighbor_i).n
      r += _epsilon * 0.5 * _phi[_j][_qp] * _mu_neighbor[_qp] *
           (_grad_test_neighbor[_i][_qp] * _normals[_qp]);

      // d/du_j of: -sigma/h * (u - u_neighbor) * v_neighbor_i
      r -= _sigma * mu_avg / h_elem * _phi[_j][_qp] * _test_neighbor[_i][_qp];
      break;

    case Moose::NeighborNeighbor:
      // d/du_neighbor_j of: +0.5 * mu_neighbor * grad(u_neighbor).n * v_neighbor_i
      r += 0.5 * _mu_neighbor[_qp] * (_grad_phi_neighbor[_j][_qp] * _normals[_qp]) *
           _test_neighbor[_i][_qp];

      // d/du_neighbor_j of: +epsilon * 0.5 * (u - u_neighbor) * mu_neighbor * grad(v_neighbor_i).n
      r += _epsilon * 0.5 * (-_phi_neighbor[_j][_qp]) * _mu_neighbor[_qp] *
           (_grad_test_neighbor[_i][_qp] * _normals[_qp]);

      // d/du_neighbor_j of: -sigma/h * (u - u_neighbor) * v_neighbor_i
      r -= _sigma * mu_avg / h_elem * (-_phi_neighbor[_j][_qp]) * _test_neighbor[_i][_qp];
      break;
  }

  return r;
}
