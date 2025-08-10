//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ElkADPoroMechanicsCoupling.h"

// MOOSE includes
#include "Function.h"
#include "MooseMesh.h"
#include "MooseVariable.h"

/**
 * Created by Chunhui Zhao, Jul 5th, 2024
 * ElkADPoroMechanicsCoupling computes -coefficient*porepressure*grad_test[component]
 */

registerMooseObject("farmsApp", ElkADPoroMechanicsCoupling);

InputParameters
ElkADPoroMechanicsCoupling::validParams()
{
  InputParameters params = ADKernel::validParams();
  params.addClassDescription(
      "Adds $-Bi \\cdot p_s \\cdot \\nabla \\Psi_c$, where the subscript $c$ is the component.");
  params.addRequiredCoupledVar("porepressure", "Pore pressure, $p_s$.");
  params.addRequiredParam<unsigned int>("component",
                                        "The gradient direction (0 for x, 1 for y and 2 for z)");
  params.addRequiredParam<bool>("multiply_biot_coefficient","whether or not multiply biot coefficient with the weak form");
  return params;
}

ElkADPoroMechanicsCoupling::ElkADPoroMechanicsCoupling(const InputParameters & parameters)
  : ADKernel(parameters),
    _multiply_biot_coefficient(getParam<bool>("multiply_biot_coefficient")),
    _coefficient(getADMaterialProperty<Real>("biot_coefficient")),
    _porepressure(adCoupledValue("porepressure")),
    _porepressure_var_num(coupled("porepressure")),
    _component(getParam<unsigned int>("component"))
{
  if (_component >= _mesh.dimension())
    mooseError("ElkADPoroMechanicsCoupling: component should not be greater than the mesh dimension\n");
}

ADReal
ElkADPoroMechanicsCoupling::computeQpResidual()
{
  if (_multiply_biot_coefficient){
    return -_coefficient[_qp] * _porepressure[_qp] * _grad_test[_i][_qp](_component);
  }
  else{
    return -_porepressure[_qp] * _grad_test[_i][_qp](_component);
  }
}