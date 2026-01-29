//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DGElasticityDirichletBC.h"
#include "Function.h"
#include "MooseVariableFE.h"
#include "libmesh/utility.h"

registerMooseObject("farmsApp", DGElasticityDirichletBC);

InputParameters
DGElasticityDirichletBC::validParams()
{
  InputParameters params = IntegratedBC::validParams();
  params.addClassDescription(
      "Dirichlet boundary condition for DG elasticity using SIPG method. "
      "Weakly enforces prescribed displacement at boundaries.");
  params.addParam<Real>("value", 0.0, "The prescribed displacement value at the boundary");
  params.addParam<FunctionName>("function", "Optional function for prescribed displacement");
  params.addParam<Real>("epsilon", 1.0,
                        "Symmetry parameter: 1 = SIPG, -1 = NIPG, 0 = IIPG");
  params.addParam<Real>("sigma", 6.0, "Penalty parameter multiplier for stability");
  params.addParam<MaterialPropertyName>("shear_modulus", "shear_modulus",
                                        "The shear modulus material property name");
  return params;
}

DGElasticityDirichletBC::DGElasticityDirichletBC(const InputParameters & parameters)
  : IntegratedBC(parameters),
    _value(getParam<Real>("value")),
    _func(isParamValid("function") ? &getFunction("function") : nullptr),
    _epsilon(getParam<Real>("epsilon")),
    _sigma(getParam<Real>("sigma")),
    _mu(getMaterialProperty<Real>("shear_modulus"))
{
}

Real
DGElasticityDirichletBC::prescribedValue() const
{
  if (_func)
    return _func->value(_t, _q_point[_qp]);
  return _value;
}

Real
DGElasticityDirichletBC::computeQpResidual()
{
  // Compute characteristic element size h
  const int elem_b_order = std::max(libMesh::Order(1), _var.order());
  const Real h_elem =
      _current_elem_volume / _current_side_volume * 1.0 / Utility::pow<2>(elem_b_order);

  const Real g = prescribedValue();
  Real r = 0.0;

  // Consistency term: -mu*grad(u).n * v
  r -= _mu[_qp] * (_grad_u[_qp] * _normals[_qp]) * _test[_i][_qp];

  // Symmetry term: epsilon * (u - g) * mu * grad(v).n
  r += _epsilon * (_u[_qp] - g) * _mu[_qp] * (_grad_test[_i][_qp] * _normals[_qp]);

  // Penalty term: sigma * mu / h * (u - g) * v
  r += _sigma * _mu[_qp] / h_elem * (_u[_qp] - g) * _test[_i][_qp];

  return r;
}

Real
DGElasticityDirichletBC::computeQpJacobian()
{
  // Compute characteristic element size h
  const int elem_b_order = std::max(libMesh::Order(1), _var.order());
  const Real h_elem =
      _current_elem_volume / _current_side_volume * 1.0 / Utility::pow<2>(elem_b_order);

  Real r = 0.0;

  // d/du of consistency term: -mu*grad(phi_j).n * v_i
  r -= _mu[_qp] * (_grad_phi[_j][_qp] * _normals[_qp]) * _test[_i][_qp];

  // d/du of symmetry term: epsilon * phi_j * mu * grad(v_i).n
  r += _epsilon * _phi[_j][_qp] * _mu[_qp] * (_grad_test[_i][_qp] * _normals[_qp]);

  // d/du of penalty term: sigma * mu / h * phi_j * v_i
  r += _sigma * _mu[_qp] / h_elem * _phi[_j][_qp] * _test[_i][_qp];

  return r;
}
