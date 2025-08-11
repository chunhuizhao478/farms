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
 * Implement flux * grad(test)
 * flux: k / mu_f grad(p)
 * k: permeability
 * mu_f: fluid viscosity
 * p: pore pressure
 */
#include "ElkADPorousFlowFullySaturatedDarcyBase.h"
#include "SubProblem.h"

registerMooseObject("farmsApp", ElkADPorousFlowFullySaturatedDarcyBase);

InputParameters ElkADPorousFlowFullySaturatedDarcyBase::validParams()
{
    InputParameters params = ADKernel::validParams();
    return params;
}

ElkADPorousFlowFullySaturatedDarcyBase::ElkADPorousFlowFullySaturatedDarcyBase(const InputParameters &
parameters):ADKernel(parameters),
    _viscosity(getADMaterialProperty<Real>("viscosity")),
    _permeability(getADMaterialProperty<RealTensorValue>("permeability"))    
{}

ADReal
ElkADPorousFlowFullySaturatedDarcyBase::computeQpResidual()
{
   return _grad_test[_i][_qp] * (_permeability[_qp] * (_grad_u[_qp])) / _viscosity[_qp];
}