//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ElasticTractionAux.h"

registerMooseObject("farmsApp", ElasticTractionAux);

InputParameters
ElasticTractionAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription(
      "Computes elastic shear traction from displacement gradient at a boundary. "
      "For antiplane shear: τ = μ * ∂w/∂n. "
      "Designed for Elasticity MainApp in MultiApp SEAS architecture.");

  params.addRequiredCoupledVar("displacement", "The displacement variable");
  params.addRequiredParam<Real>("shear_modulus", "Shear modulus μ (Pa)");
  params.addParam<unsigned int>(
      "normal_component", 0, "Normal direction component (0=x, 1=y, 2=z), default x");

  return params;
}

ElasticTractionAux::ElasticTractionAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _disp(coupledValue("displacement")),
    _grad_disp(coupledGradient("displacement")),
    _shear_modulus(getParam<Real>("shear_modulus")),
    _normal_component(getParam<unsigned int>("normal_component"))
{
  if (_shear_modulus <= 0.0)
    mooseError("Shear modulus must be positive");
  if (_normal_component > 2)
    mooseError("Normal component must be 0 (x), 1 (y), or 2 (z)");
}

Real
ElasticTractionAux::computeValue()
{
  // For antiplane shear: τ = μ * ∂w/∂n
  // The gradient component in the normal direction gives the shear strain
  Real grad_n = _grad_disp[_qp](_normal_component);

  return _shear_modulus * grad_n;
}
