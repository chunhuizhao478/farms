//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ElkPorousFlowFullySaturatedMassTimeDerivative.h"

#include "MooseVariable.h"

registerMooseObject("farmsApp", ElkPorousFlowFullySaturatedMassTimeDerivative);

InputParameters
ElkPorousFlowFullySaturatedMassTimeDerivative::validParams()
{
  InputParameters params = TimeKernel::validParams();
  MooseEnum coupling_type("Hydro ThermoHydro HydroMechanical ThermoHydroMechanical", "Hydro");
  params.addParam<MooseEnum>(
      "coupling_type",
      coupling_type,
      "The type of simulation. For Mechanical, supply Biot coefficient. For Thermal, supply thermal expansion material");
  params.addRangeCheckedParam<Real>(
      "biot_coefficient", 1.0, "biot_coefficient>=0 & biot_coefficient<=1", "Biot coefficient (constant, ignored if use_damaged_biot=true)");
  params.addParam<bool>("use_damaged_biot", false,
                        "Use biot_coefficient from material property 'biot_coefficient'");
  params.addParam<bool>("multiply_by_density",
                        true,
                        "If true, then this Kernel is the time derivative of the fluid mass; otherwise, derivative of the fluid volume");
  params.addRequiredParam<UserObjectName>(
      "PorousFlowDictator", "The UserObject that holds the list of PorousFlow variable names.");
  params.addClassDescription(
      "Fully-saturated single-component, single-phase fluid mass derivative wrt time (with optional damaged Biot coefficient)");
  return params;
}

ElkPorousFlowFullySaturatedMassTimeDerivative::ElkPorousFlowFullySaturatedMassTimeDerivative(
    const InputParameters & parameters)
  : TimeKernel(parameters),
    _dictator(getUserObject<PorousFlowDictator>("PorousFlowDictator")),
    _var_is_porflow_var(_dictator.isPorousFlowVariable(_var.number())),
    _multiply_by_density(getParam<bool>("multiply_by_density")),
    _coupling_type(getParam<MooseEnum>("coupling_type").getEnum<CouplingTypeEnum>()),
    _includes_thermal(_coupling_type == CouplingTypeEnum::ThermoHydro ||
                      _coupling_type == CouplingTypeEnum::ThermoHydroMechanical),
    _includes_mechanical(_coupling_type == CouplingTypeEnum::HydroMechanical ||
                         _coupling_type == CouplingTypeEnum::ThermoHydroMechanical),
    _biot_coefficient_const(getParam<Real>("biot_coefficient")),
    _use_damaged_biot(getParam<bool>("use_damaged_biot")),
    _biot_coefficient_mp(_use_damaged_biot ? &getMaterialProperty<Real>("biot_coefficient")
                                           : nullptr),
    _biot_modulus(getMaterialProperty<Real>("PorousFlow_constant_biot_modulus_qp")),
    _thermal_coeff(_includes_thermal ? &getMaterialProperty<Real>(
                                       "PorousFlow_constant_thermal_expansion_coefficient_qp")
                                     : nullptr),
    _fluid_density(_multiply_by_density ? &getMaterialProperty<std::vector<Real>>(
                                              "PorousFlow_fluid_phase_density_qp")
                                        : nullptr),
    _dfluid_density_dvar(_multiply_by_density
                             ? &getMaterialProperty<std::vector<std::vector<Real>>>(
                                   "dPorousFlow_fluid_phase_density_qp_dvar")
                             : nullptr),
    _pp(getMaterialProperty<std::vector<Real>>("PorousFlow_porepressure_qp")),
    _pp_old(getMaterialPropertyOld<std::vector<Real>>("PorousFlow_porepressure_qp")),
    _dpp_dvar(
        getMaterialProperty<std::vector<std::vector<Real>>>("dPorousFlow_porepressure_qp_dvar")),
    _temperature(_includes_thermal ? &getMaterialProperty<Real>("PorousFlow_temperature_qp")
                                   : nullptr),
    _temperature_old(_includes_thermal ? &getMaterialPropertyOld<Real>("PorousFlow_temperature_qp")
                                       : nullptr),
    _dtemperature_dvar(_includes_thermal
                           ? &getMaterialProperty<std::vector<Real>>(
                                 "dPorousFlow_temperature_qp_dvar")
                           : nullptr),
    _strain_rate(_includes_mechanical ? &getMaterialProperty<Real>(
                                         "PorousFlow_volumetric_strain_rate_qp")
                                      : nullptr),
    _dstrain_rate_dvar(_includes_mechanical ? &getMaterialProperty<std::vector<RealGradient>>(
                                                  "dPorousFlow_volumetric_strain_rate_qp_dvar")
                                            : nullptr)
{
  if (_dictator.numComponents() != 1 || _dictator.numPhases() != 1)
    mooseError("ElkPorousFlowFullySaturatedMassTimeDerivative is only applicable to single-phase, "
               "single-component fluid-flow problems.");
}

Real
ElkPorousFlowFullySaturatedMassTimeDerivative::computeQpResidual()
{
  const unsigned phase = 0;
  const Real alpha = biot();
  Real volume = (_pp[_qp][phase] - _pp_old[_qp][phase]) / _dt / _biot_modulus[_qp];
  if (_includes_thermal)
    volume -= (*_thermal_coeff)[_qp] * ((*_temperature)[_qp] - (*_temperature_old)[_qp]) / _dt;
  if (_includes_mechanical)
    volume += alpha * (*_strain_rate)[_qp];
  if (_multiply_by_density)
    return _test[_i][_qp] * (*_fluid_density)[_qp][phase] * volume;
  return _test[_i][_qp] * volume;
}

Real
ElkPorousFlowFullySaturatedMassTimeDerivative::computeQpJacobian()
{
  // If the variable is not a PorousFlow variable (very unusual), the diag Jacobian terms are 0
  if (!_var_is_porflow_var)
    return 0.0;
  return computeQpJac(_dictator.porousFlowVariableNum(_var.number()));
}

Real
ElkPorousFlowFullySaturatedMassTimeDerivative::computeQpOffDiagJacobian(unsigned int jvar)
{
  // If the variable is not a PorousFlow variable, the OffDiag Jacobian terms are 0
  if (_dictator.notPorousFlowVariable(jvar))
    return 0.0;
  return computeQpJac(_dictator.porousFlowVariableNum(jvar));
}

Real
ElkPorousFlowFullySaturatedMassTimeDerivative::computeQpJac(unsigned int pvar)
{
  const unsigned phase = 0;
  const Real alpha = biot();
  Real volume = (_pp[_qp][phase] - _pp_old[_qp][phase]) / _dt / _biot_modulus[_qp];
  Real dvolume = _dpp_dvar[_qp][phase][pvar] / _dt / _biot_modulus[_qp] * _phi[_j][_qp];
  if (_includes_thermal)
  {
    volume -= (*_thermal_coeff)[_qp] * ((*_temperature)[_qp] - (*_temperature_old)[_qp]) / _dt;
    dvolume -= (*_thermal_coeff)[_qp] * (*_dtemperature_dvar)[_qp][pvar] / _dt * _phi[_j][_qp];
  }
  if (_includes_mechanical)
  {
    volume += alpha * (*_strain_rate)[_qp];
    dvolume += alpha * (*_dstrain_rate_dvar)[_qp][pvar] * _grad_phi[_j][_qp];
  }
  if (_multiply_by_density)
    return _test[_i][_qp] * ((*_fluid_density)[_qp][phase] * dvolume +
                             (*_dfluid_density_dvar)[_qp][phase][pvar] * _phi[_j][_qp] * volume);
  return _test[_i][_qp] * dvolume;
}
