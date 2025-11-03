//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ElkPorousFlowEffectiveStressCoupling.h"

#include "Function.h"
#include "MooseMesh.h"
#include "MooseVariable.h"

registerMooseObject("farmsApp", ElkPorousFlowEffectiveStressCoupling);

InputParameters
ElkPorousFlowEffectiveStressCoupling::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Implements the weak form of alpha * grad(effective fluid pressure) with optional damaged Biot coefficient");
  params.addRequiredParam<UserObjectName>(
      "PorousFlowDictator", "The UserObject that holds the list of PorousFlow variable names.");
  params.addRangeCheckedParam<Real>(
      "biot_coefficient", 1, "biot_coefficient>=0&biot_coefficient<=1", "Biot coefficient (constant, ignored if use_damaged_biot=true)");
  params.addRequiredParam<unsigned int>("component",
                                        "The component (0 for x, 1 for y and 2 for z) of grad(P)");
  params.addParam<bool>("use_damaged_biot", false, "Use biot_coefficient from material property 'biot_coefficient_damaged'");
  return params;
}

ElkPorousFlowEffectiveStressCoupling::ElkPorousFlowEffectiveStressCoupling(
    const InputParameters & parameters)
  : Kernel(parameters),
    _dictator(getUserObject<PorousFlowDictator>("PorousFlowDictator")),
    _coefficient_const(getParam<Real>("biot_coefficient")),
    _use_damaged_biot(getParam<bool>("use_damaged_biot")),
    _biot_coeff_mp(_use_damaged_biot ? &getMaterialProperty<Real>("biot_coefficient_damaged") : nullptr),
    _component(getParam<unsigned int>("component")),
    _pf(getMaterialProperty<Real>("PorousFlow_effective_fluid_pressure_qp")),
    _dpf_dvar(
        getMaterialProperty<std::vector<Real>>("dPorousFlow_effective_fluid_pressure_qp_dvar")),
    _rz(getBlockCoordSystem() == Moose::COORD_RZ)
{
  if (_component >= _mesh.dimension())
    paramError("component", "The component cannot be greater than the mesh dimension");
}

Real
ElkPorousFlowEffectiveStressCoupling::computeQpResidual()
{
  const Real alpha = biot();
  if (_rz && _component == 0)
    return -alpha * _pf[_qp] * (_grad_test[_i][_qp](0) + _test[_i][_qp] / _q_point[_qp](0));
  return -alpha * _pf[_qp] * _grad_test[_i][_qp](_component);
}

Real
ElkPorousFlowEffectiveStressCoupling::computeQpJacobian()
{
  const Real alpha = biot();
  if (_dictator.notPorousFlowVariable(_var.number()))
    return 0.0;
  const unsigned int pvar = _dictator.porousFlowVariableNum(_var.number());
  if (_rz && _component == 0)
    return -alpha * _phi[_j][_qp] * _dpf_dvar[_qp][pvar] *
           (_grad_test[_i][_qp](0) + _test[_i][_qp] / _q_point[_qp](0));
  return -alpha * _phi[_j][_qp] * _dpf_dvar[_qp][pvar] * _grad_test[_i][_qp](_component);
}

Real
ElkPorousFlowEffectiveStressCoupling::computeQpOffDiagJacobian(unsigned int jvar)
{
  const Real alpha = biot();
  if (_dictator.notPorousFlowVariable(jvar))
    return 0.0;
  const unsigned int pvar = _dictator.porousFlowVariableNum(jvar);
  if (_rz && _component == 0)
    return -alpha * _phi[_j][_qp] * _dpf_dvar[_qp][pvar] *
           (_grad_test[_i][_qp](0) + _test[_i][_qp] / _q_point[_qp](0));
  return -alpha * _phi[_j][_qp] * _dpf_dvar[_qp][pvar] * _grad_test[_i][_qp](_component);
}
