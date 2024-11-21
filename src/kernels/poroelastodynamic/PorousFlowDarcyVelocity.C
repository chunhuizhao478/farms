#include "PorousFlowDarcyVelocity.h"
#include "MooseVariable.h"


registerMooseObject("farmsApp", PorousFlowDarcyVelocity);

InputParameters
PorousFlowDarcyVelocity::validParams()
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
  params.addRequiredRangeCheckedParam<unsigned int>("component",
                                                    "component < 3",
                                                    "An integer corresponding to the direction "
                                                    "the variable this kernel acts in. (0 for x, "
                                                    "1 for y, 2 for z)");
  return params;
}

PorousFlowDarcyVelocity::PorousFlowDarcyVelocity(
    const InputParameters & parameters)
  : Kernel(parameters),
    _multiply_by_density(getParam<bool>("multiply_by_density")),
    _permeability(getMaterialProperty<RealTensorValue>("PorousFlow_permeability_qp")),
    _dpermeability_dvar(
        getMaterialProperty<std::vector<RealTensorValue>>("dPorousFlow_permeability_qp_dvar")),
    _dpermeability_dgradvar(getMaterialProperty<std::vector<std::vector<RealTensorValue>>>(
        "dPorousFlow_permeability_qp_dgradvar")),
    _density(getMaterialProperty<std::vector<Real>>("PorousFlow_fluid_phase_density_qp")),
    _ddensity_dvar(getMaterialProperty<std::vector<std::vector<Real>>>(
        "dPorousFlow_fluid_phase_density_qp_dvar")),
    _viscosity(getMaterialProperty<std::vector<Real>>("PorousFlow_viscosity_qp")),
    _dviscosity_dvar(
        getMaterialProperty<std::vector<std::vector<Real>>>("dPorousFlow_viscosity_qp_dvar")),
    _dictator(getUserObject<PorousFlowDictator>("PorousFlowDictator")),
    _component(getParam<unsigned int>("component")),
    _perm_derivs(_dictator.usePermDerivs())
{
  if (_dictator.numPhases() != 1)
    mooseError("PorousFlowFullySaturatedDarcyBase should not be used for multi-phase scenarios as "
               "it does no upwinding and does not include relative-permeability effects");
}

Real
PorousFlowDarcyVelocity::computeQpResidual()
{
  const unsigned ph = 0;
  const Real mob = mobility();
      
  return _test[_i][_qp]  * mob * _u[_qp] / _permeability[_qp](_component,_component);

}

Real
PorousFlowDarcyVelocity::computeQpJacobian()
{
  return computeQpOffDiagJacobian(_var.number());
}

Real
PorousFlowDarcyVelocity::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_dictator.notPorousFlowVariable(jvar))
    return 0.0;

  const unsigned ph = 0;
  const unsigned pvar = _dictator.porousFlowVariableNum(jvar);
  const Real mob = mobility();
  const Real dmob = dmobility(pvar) ;
  
   const RealVectorValue flow = _u[_qp] / _permeability[_qp](_component,_component);
    
  RealVectorValue dflow = _phi[_j][_qp] / _permeability[_qp](_component,_component);


 // if (_perm_derivs)
 // {
 //   dflow -=   _u[_qp] * _dpermeability_dvar[_qp][pvar](_component,_component) * _phi[_j][_qp] 
 //             /std::pow(_permeability[_qp](_component,_component),2);
             
  //  for (const auto i : make_range(Moose::dim))
  //    dflow -=  _u[_qp] * _dpermeability_dgradvar[_qp][i][pvar](_component,_component)  * _grad_phi[_j][_qp](i) 
   //               /std::pow(_permeability[_qp](_component,_component),2);
  //}
  
  return _test[_i][_qp]  * mob * _phi[_j][_qp] / _permeability[_qp](_component,_component) + 
         _test[_i][_qp]  *  _u[_qp] / _permeability[_qp](_component,_component) * dmob ;
}
Real
PorousFlowDarcyVelocity::mobility() const
{
  const unsigned ph = 0;
  Real mob =  _viscosity[_qp][ph]  ;

  if (_multiply_by_density)
    mob *= _density[_qp][ph];
  return mob;
}

Real
PorousFlowDarcyVelocity::dmobility(unsigned pvar) const
{
  const unsigned ph = 0;
  Real mob =  _viscosity[_qp][ph];
  Real dmob = _dviscosity_dvar[_qp][ph][pvar];


  if (_multiply_by_density)
    dmob = _density[_qp][ph] * dmob + _ddensity_dvar[_qp][ph][pvar] * mob ;

  return dmob;
}
 