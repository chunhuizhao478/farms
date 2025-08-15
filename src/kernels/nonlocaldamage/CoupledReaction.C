//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CoupledReaction.h"

registerMooseObject("farmsApp", CoupledReaction);

InputParameters
CoupledReaction::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription(
      "Implements a simple consuming reaction term with weak form $(\\psi_i, \\lambda u_h)$.");
  params.addParam<Real>(
      "rate", 1.0, "The $(\\lambda)$ multiplier, the relative amount consumed per unit time.");
  params.declareControllable("rate");
  params.addRequiredCoupledVar(
      "eqstrain_local",
      "The local equivalent strain used in the damage evolution law");
  params.addRequiredParam<Real>(
      "length_scale",
      "The length scale used in the gradient activity parameter for the equivalent strain");
  params.addRequiredParam<Real>(
      "kappa_i",
      "The equivalent strain at which the gradient activity starts");
  params.addRequiredParam<Real>(
      "c0",
      "The minimum value of the gradient activity parameter for the equivalent strain");
  return params;
}

CoupledReaction::CoupledReaction(const InputParameters & parameters)
  : Kernel(parameters), 
    _rate(getParam<Real>("rate")), 
    _eqstrain_local(coupledValue("eqstrain_local")),
    _length_scale(getParam<Real>("length_scale")),
    _kappa_i(getParam<Real>("kappa_i")),
    _c0(getParam<Real>("c0"))
{
}

Real
CoupledReaction::computeQpResidual()
{
  // match residual scaling: (_rate/xi) * u
  Real xi = 0.0;
  const Real l = _length_scale;
  const Real c = 0.2 * l * l;
  const Real e_xi = 10.0 * _kappa_i;
  const Real e = std::max(_eqstrain_local[_qp], 0.0);
  xi = (e < e_xi) ? _c0 + (c - _c0) * (e / e_xi) : c;

  return _test[_i][_qp] * _rate * _u[_qp] / xi;
}

Real
CoupledReaction::computeQpJacobian()
{

  // match residual scaling: (_rate/xi) * u
  Real xi = 0.0;
  const Real l = _length_scale;
  const Real c = 0.2 * l * l;
  const Real e_xi = 10.0 * _kappa_i;
  const Real e = std::max(_eqstrain_local[_qp], 0.0);
  xi = (e < e_xi) ? _c0 + (c - _c0) * (e / e_xi) : c;

  return _test[_i][_qp] * _rate / xi * _phi[_j][_qp];
}
