#include "PorousFlowFullySaturatedMassTimeDerivativeUserHM.h"

#include "MooseVariable.h"

registerMooseObject("farmsApp", PorousFlowFullySaturatedMassTimeDerivativeUserHM);

InputParameters
PorousFlowFullySaturatedMassTimeDerivativeUserHM::validParams()
{
  InputParameters params = TimeKernel::validParams();
  MooseEnum coupling_type("HydroMechanical", "Hydro");
  params.addParam<MooseEnum>("coupling_type",
                             coupling_type,
                             "The type of simulation.  For simulations involving Mechanical "
                             "deformations");
  params.addRangeCheckedParam<Real>(
      "biot_coefficient", 1.0, "biot_coefficient>=0 & biot_coefficient<=1", "Biot coefficient");
  params.addParam<bool>("multiply_by_density",
                        true,
                        "If true, then this Kernel is the time derivative of the fluid "
                        "mass.  If false, then this Kernel is the derivative of the "
                        "fluid volume (which is common in poro-mechanics)");
  params.addRequiredParam<UserObjectName>(
      "PorousFlowDictator", "The UserObject that holds the list of PorousFlow variable names.");
  params.addClassDescription("Fully-saturated version of the single-component, single-phase fluid "
                             "mass derivative wrt time");
  return params;
}

PorousFlowFullySaturatedMassTimeDerivativeUserHM::PorousFlowFullySaturatedMassTimeDerivativeUserHM(
    const InputParameters & parameters)
  : TimeKernel(parameters),
    _dictator(getUserObject<PorousFlowDictator>("PorousFlowDictator")),
    _var_is_porflow_var(_dictator.isPorousFlowVariable(_var.number())),
    _multiply_by_density(getParam<bool>("multiply_by_density")),
    _coupling_type(getParam<MooseEnum>("coupling_type").getEnum<CouplingTypeEnum>()),
    _includes_mechanical(_coupling_type == CouplingTypeEnum::HydroMechanical),
    _biot_coefficient(getParam<Real>("biot_coefficient")),
    _biot_modulus(getMaterialProperty<Real>("PorousFlow_constant_biot_modulus_qp")),
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
{
  if (_dictator.numComponents() != 1 || _dictator.numPhases() != 1)
    mooseError("PorousFlowFullySaturatedMassTimeDerivativeUserHM is only applicable to single-phase, "
               "single-component fluid-flow problems.  The Dictator proclaims that you have more "
               "than one phase or more than one fluid component.  The Dictator does not take such "
               "mistakes lightly");
}

Real
PorousFlowFullySaturatedMassTimeDerivativeUserHM::computeQpResidual()
{
  const unsigned phase = 0;
  Real volume = (_pp[_qp][phase] - _pp_old[_qp][phase]) / _dt / _biot_modulus[_qp]/3;
  if (_includes_mechanical)
    volume += _biot_coefficient * (*_u_dot)[_qp]*1/3;
  if (_multiply_by_density)
    return _test[_i][_qp] * (*_fluid_density)[_qp][phase] * volume;
  return _test[_i][_qp] * volume;
}

Real
PorousFlowFullySaturatedMassTimeDerivativeUserHM::computeQpJacobian()
{
  // If the variable is not a PorousFlow variable (very unusual), the diag Jacobian terms are 0
  if (!_var_is_porflow_var)
    return 0.0;
  return computeQpJac(_dictator.porousFlowVariableNum(_var.number()));
}

Real
PorousFlowFullySaturatedMassTimeDerivativeUserHM::computeQpOffDiagJacobian(unsigned int jvar)
{
  // If the variable is not a PorousFlow variable, the OffDiag Jacobian terms are 0
  if (_dictator.notPorousFlowVariable(jvar))
    return 0.0;
  return computeQpJac(_dictator.porousFlowVariableNum(jvar));
}

Real
PorousFlowFullySaturatedMassTimeDerivativeUserHM::computeQpJac(unsigned int pvar)
{
  const unsigned phase = 0;
  Real volume = (_pp[_qp][phase] - _pp_old[_qp][phase]) / _dt / _biot_modulus[_qp]/3;
  Real dvolume = _dpp_dvar[_qp][phase][pvar] / _dt / _biot_modulus[_qp] * _phi[_j][_qp]/3;
  if (_includes_mechanical)
  {
    volume += _biot_coefficient * (*_u_dot)[_qp]*1/3;
    dvolume += _biot_coefficient * (*_du_dot_du)[_qp][pvar] * _grad_phi[_j][_qp]*1/3;
  }
  if (_multiply_by_density)
    return _test[_i][_qp] * ((*_fluid_density)[_qp][phase] * dvolume +
                             (*_dfluid_density_dvar)[_qp][phase][pvar] * _phi[_j][_qp] * volume);
  return _test[_i][_qp] * dvolume;
}