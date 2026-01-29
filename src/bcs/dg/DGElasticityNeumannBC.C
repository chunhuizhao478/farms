//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DGElasticityNeumannBC.h"
#include "Function.h"

registerMooseObject("farmsApp", DGElasticityNeumannBC);

InputParameters
DGElasticityNeumannBC::validParams()
{
  InputParameters params = IntegratedBC::validParams();
  params.addClassDescription(
      "Neumann (traction) boundary condition for DG elasticity. "
      "Applies prescribed traction at boundaries. Use traction=0 for free surface.");
  params.addParam<Real>("traction", 0.0, "The prescribed traction value at the boundary");
  params.addParam<FunctionName>("function", "Optional function for prescribed traction");
  return params;
}

DGElasticityNeumannBC::DGElasticityNeumannBC(const InputParameters & parameters)
  : IntegratedBC(parameters),
    _traction(getParam<Real>("traction")),
    _func(isParamValid("function") ? &getFunction("function") : nullptr)
{
}

Real
DGElasticityNeumannBC::computeQpResidual()
{
  // Get prescribed traction value
  Real t = _traction;
  if (_func)
    t = _func->value(_t, _q_point[_qp]);

  // Residual: -t * v (traction applied to boundary)
  // This comes from the weak form: ∫_Γ σ·n * v dA = ∫_Γ t * v dA
  // Moving to RHS gives -t * v
  return -t * _test[_i][_qp];
}

Real
DGElasticityNeumannBC::computeQpJacobian()
{
  // Neumann BC has no dependence on the solution
  return 0.0;
}
