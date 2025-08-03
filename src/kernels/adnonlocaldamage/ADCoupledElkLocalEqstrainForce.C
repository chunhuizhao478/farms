//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADCoupledElkLocalEqstrainForce.h"

registerADMooseObject("farmsApp", ADCoupledElkLocalEqstrainForce);

InputParameters
ADCoupledElkLocalEqstrainForce::validParams()
{
  InputParameters params = ADKernel::validParams();
  params.addClassDescription("AD Kernel for implement local equivalent strain force");
  // params.addRequiredCoupledVar(
  //     "eqstrain_local",
  //     "The local equivalent strain used in the damage evolution law");
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

ADCoupledElkLocalEqstrainForce::ADCoupledElkLocalEqstrainForce(const InputParameters & parameters)
  : ADKernel(parameters),
    _eqstrain_local(getADMaterialProperty<Real>("eqstrain_local")),
    _length_scale(getParam<Real>("length_scale")),
    _kappa_i(getParam<Real>("kappa_i")),
    _c0(getParam<Real>("c0"))
{
}

ADReal
ADCoupledElkLocalEqstrainForce::computeQpResidual()
{
  // Test gradient activity parameter
  ADReal xi = 0.0;
  Real l = _length_scale; // Length scale
  Real kappa_i = _kappa_i;
  Real e_xi = 10*kappa_i; // xi is the equivalent strain at which the gradient activity starts
  Real c0 = _c0;
  Real c = 0.5*l*l;
  
  if (_eqstrain_local[_qp] < e_xi){
    xi = c0 + (c - c0) * (_eqstrain_local[_qp] / e_xi);
  }
  else{
    xi = c;
  }

  return -1.0 * _test[_i][_qp] * _eqstrain_local[_qp] / xi;
}