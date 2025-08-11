//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

/**
 * Created by Chunhui Zhao, Aug 11th, 2025
 * Implement alpha * volumetric strain rate * test + 1 / biot modulus * pressure rate * test
 * alpha: biot coefficient
 */
#include "ElkADPorousFlowFullySaturatedMassTimeDerivative.h"
#include "SubProblem.h"

registerMooseObject("farmsApp", ElkADPorousFlowFullySaturatedMassTimeDerivative);

InputParameters ElkADPorousFlowFullySaturatedMassTimeDerivative::validParams()
{
    InputParameters params = ADTimeKernelValue::validParams();
    return params;
}

ElkADPorousFlowFullySaturatedMassTimeDerivative::ElkADPorousFlowFullySaturatedMassTimeDerivative(const InputParameters &
parameters)
    : ADTimeKernelValue(parameters),
        _biot_modulus(getADMaterialProperty<Real>("biot_modulus")),
        _biot_coefficient(getADMaterialProperty<Real>("biot_coefficient")),
        _vol_strain(getADMaterialProperty<Real>("vol_strain")),
        _vol_strain_old(getMaterialPropertyOldByName<Real>("vol_strain"))
{}

ADReal
ElkADPorousFlowFullySaturatedMassTimeDerivative::precomputeQpResidual()
{
    // Pressure rate term: p_dot / M  (M = biot_modulus)
    ADReal r = _u_dot[_qp] / _biot_modulus[_qp];
    // Volumetric strain rate contribution: alpha * d(eps_v)/dt
    r += _biot_coefficient[_qp] * (_vol_strain[_qp] - _vol_strain_old[_qp]) / _dt;
    return r;
}