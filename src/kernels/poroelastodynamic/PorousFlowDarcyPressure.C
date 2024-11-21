

#include "PorousFlowDarcyPressure.h"

#include "MooseVariable.h"

registerMooseObject("farmsApp", PorousFlowDarcyPressure);

InputParameters
PorousFlowDarcyPressure::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addParam<bool>("multiply_by_density",
                        false,
                        "If true, then this Kernel is the fluid mass "
                        "flux.  If false, then this Kernel is the "
                        "fluid volume flux (which is common in "
                        "poro-mechanics)");
  params.addRequiredParam<UserObjectName>(
      "PorousFlowDictator", "The UserObject that holds the list of PorousFlow variable names");
  params.addClassDescription("Darcy flux suitable for models involving a fully-saturated, single "
                             "phase, single component fluid.  No upwinding is used");
  params.addRequiredParam<unsigned int>("component",
                                        "The component (0 for x, 1 for y and 2 for z) of P");
  return params;
}

PorousFlowDarcyPressure::PorousFlowDarcyPressure(
    const InputParameters & parameters)
  : Kernel(parameters),
    _multiply_by_density(getParam<bool>("multiply_by_density")),
    _density(getMaterialProperty<std::vector<Real>>("PorousFlow_fluid_phase_density_qp")),
    _ddensity_dvar(getMaterialProperty<std::vector<std::vector<Real>>>(
        "dPorousFlow_fluid_phase_density_qp_dvar")),
    _pp(getMaterialProperty<std::vector<Real>>("PorousFlow_porepressure_qp")),
    _dpp_dvar(getMaterialProperty<std::vector<std::vector<Real>>>("dPorousFlow_porepressure_qp_dvar")),
    _dictator(getUserObject<PorousFlowDictator>("PorousFlowDictator")),
    _component(getParam<unsigned int>("component"))
{
}

Real
PorousFlowDarcyPressure::computeQpResidual()
{
  const unsigned ph = 0;

  if (_multiply_by_density)
    return   _density[_qp][ph] * (- _grad_test[_i][_qp](_component) * _pp[_qp][ph]);

  return  - _grad_test[_i][_qp](_component) * _pp[_qp][ph] ;
}

Real
PorousFlowDarcyPressure::computeQpJacobian()
{
  return computeQpOffDiagJacobian(_var.number());
}

Real
PorousFlowDarcyPressure::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_dictator.notPorousFlowVariable(jvar))
    return 0.0;

  const unsigned ph = 0;
  const unsigned pvar = _dictator.porousFlowVariableNum(jvar);

  Real flow = - _grad_test[_i][_qp](_component) * _pp[_qp][ph]  ;
  Real dflow =  - _grad_test[_i][_qp](_component)* _phi[_j][_qp] * _dpp_dvar[_qp][ph][pvar];
                

  if (_multiply_by_density)
    return _density[_qp][ph] * dflow +
           _ddensity_dvar[_qp][ph][pvar] * _phi[_j][_qp] * flow;

  return  dflow;

}

