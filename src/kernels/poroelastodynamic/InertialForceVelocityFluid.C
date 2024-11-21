
#include "InertialForceVelocityFluid.h"
#include "SubProblem.h"
#include "TimeIntegrator.h"
#include "NonlinearSystemBase.h"

registerMooseObject("farmsApp", InertialForceVelocityFluid);

InputParameters
InertialForceVelocityFluid::validParams()
{
  InputParameters params = TimeKernel ::validParams();
  params.addClassDescription("Calculates the residual for the inertial force "
                             "($M \\cdot acceleration$) and the contribution of mass"
                             " dependent Rayleigh damping and HHT time "
                             " integration scheme ($\\eta \\cdot M \\cdot"
                             " ((1+\\alpha)velq2-\\alpha \\cdot vel-old) $)");
  params.set<bool>("use_displaced_mesh") = false;
  params.addParam<MaterialPropertyName>(
      "density", "density", "Name of Material Property that provides the density");
  params.addParam<unsigned int>(
      "fluid_component", 0, "The index corresponding to the component for this kernel"); 
  params.addParam<Real>(
      "apparent_density", 0); 
  params.addRequiredParam<UserObjectName>(
      "PorousFlowDictator", "The UserObject that holds the list of PorousFlow variable names.");
  return params;
}

InertialForceVelocityFluid::InertialForceVelocityFluid(const InputParameters & parameters)
  : TimeKernel(parameters),
    _fluid_component(getParam<unsigned int>("fluid_component")),
    _apparent_density(getParam<Real>("apparent_density")),
    _dictator(getUserObject<PorousFlowDictator>("PorousFlowDictator")),
    _num_phases(_dictator.numPhases()),
    _fluid_density_qp(getMaterialProperty<std::vector<Real>>("PorousFlow_fluid_phase_density_qp")),
    _porosity(getMaterialProperty<Real>("PorousFlow_porosity_qp")),
    _time_integrator(*_sys.getTimeIntegrator())
{
    _du_dot_du = &(_var.duDotDu());

    this->addFEVariableCoupleableVectorTag(_time_integrator.uDotFactorTag());
    this->addFEVariableCoupleableVectorTag(_time_integrator.uDotDotFactorTag());

    _u_dot_factor = &_var.vectorTagValue(_time_integrator.uDotFactorTag());

}

Real
InertialForceVelocityFluid::computeQpResidual()
{

  Real dens = 0.0;
  for (unsigned ph = 0; ph < _num_phases; ++ph)
  {
     dens +=  _fluid_density_qp[_qp][ph] /_porosity[_qp] + _apparent_density /_porosity[_qp] /_porosity[_qp] ;
  }

  if (_dt == 0)
    return 0;

  // Consistent mass option
  // Same for explicit, implicit, and implicit with HHT
  else
    return _test[_i][_qp]  * dens  *  (*_u_dot_factor)[_qp] ;
}


Real
InertialForceVelocityFluid::computeQpJacobian()
{
  Real dens = 0.0;
  for (unsigned ph = 0; ph < _num_phases; ++ph)
  {
     dens +=  _fluid_density_qp[_qp][ph] /_porosity[_qp] + _apparent_density /_porosity[_qp] /_porosity[_qp] ;
  }

  if (_dt == 0)
    return 0;
  else

  return _test[_i][_qp] * dens * (*_du_dot_du)[_qp] * _phi[this->_j][_qp];
}
