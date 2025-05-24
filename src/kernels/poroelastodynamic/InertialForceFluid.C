
#include "InertialForceFluid.h"
#include "SubProblem.h"
#include "TimeIntegrator.h"
#include "NonlinearSystemBase.h"

registerMooseObject("farmsApp", InertialForceFluid);

InputParameters
InertialForceFluid::validParams()
{
  InputParameters params = TimeKernel ::validParams();
  params.addClassDescription("Calculates the residual for the inertial force "
                             "($M \\cdot acceleration$) and the contribution of mass"
                             " dependent Rayleigh damping and HHT time "
                             " integration scheme ($\\eta \\cdot M \\cdot"
                             " ((1+\\alpha)velq2-\\alpha \\cdot vel-old) $)");
  params.set<bool>("use_displaced_mesh") = true;
  params.addCoupledVar("velocity", "velocity variable");
  params.addCoupledVar("acceleration", "acceleration variable");
  params.addParam<Real>("beta", "beta parameter for Newmark Time integration");
  params.addParam<Real>("gamma", "gamma parameter for Newmark Time integration");
  params.addParam<Real>("eta",
                                        0.0,
                                        "Name of material property or a constant real "
                                        "number defining the eta parameter for the "
                                        "Rayleigh damping.");
  params.addParam<MaterialPropertyName>(
      "density_scaling",
      0.0,
      "Name of material property to add mass scaling in explicit simulations.");
  params.addParam<Real>("alpha",
                        0,
                        "alpha parameter for mass dependent numerical damping induced "
                        "by HHT time integration scheme");
  params.addParam<MaterialPropertyName>(
      "density", "density", "Name of Material Property that provides the density");
  params.addParam<unsigned int>(
      "fluid_component", 0, "The index corresponding to the component for this kernel"); 
  params.addRequiredParam<UserObjectName>(
      "PorousFlowDictator", "The UserObject that holds the list of PorousFlow variable names.");
  return params;
}

InertialForceFluid::InertialForceFluid(const InputParameters & parameters)
  : TimeKernel(parameters),
    _has_beta(this->isParamValid("beta")),
    _has_gamma(this->isParamValid("gamma")),
    _beta(_has_beta ? this->template getParam<Real>("beta") : 0.1),
    _gamma(_has_gamma ? this->template getParam<Real>("gamma") : 0.1),
    _has_velocity(this->isParamValid("velocity")),
    _has_acceleration(this->isParamValid("acceleration")),
    _eta(getParam<Real>("eta")),
    _alpha(this->template getParam<Real>("alpha")),
    _fluid_component(getParam<unsigned int>("fluid_component")),
    _dictator(getUserObject<PorousFlowDictator>("PorousFlowDictator")),
    _num_phases(_dictator.numPhases()),
    _fluid_density_qp(getMaterialProperty<std::vector<Real>>("PorousFlow_fluid_phase_density_qp")),
    _porosity(getMaterialProperty<Real>("PorousFlow_porosity_qp")),
    _time_integrator(_sys.getTimeIntegrator())
{
  if (_has_beta && _has_gamma && _has_velocity && _has_acceleration)
  {
    _vel_old = &this->coupledValueOld("velocity");
    _accel_old = &this->coupledValueOld("acceleration");
    _u_old = &this->valueOld();
  }
  else if (!_has_beta && !_has_gamma && !_has_velocity && !_has_acceleration)
  {
    _u_dot_old = &(_var.uDotOld());
    _du_dot_du = &(_var.duDotDu());
    _du_dotdot_du = &(_var.duDotDotDu());

    this->addFEVariableCoupleableVectorTag(_time_integrator.uDotFactorTag());
    this->addFEVariableCoupleableVectorTag(_time_integrator.uDotDotFactorTag());

    _u_dot_factor = &_var.vectorTagValue(_time_integrator.uDotFactorTag());
    _u_dotdot_factor = &_var.vectorTagValue(_time_integrator.uDotDotFactorTag());

    if (_time_integrator.isLumped())
    {
      _u_dot_factor_dof = &_var.vectorTagDofValue(_time_integrator.uDotFactorTag());
      _u_dotdot_factor_dof = &_var.vectorTagDofValue(_time_integrator.uDotDotFactorTag());
    }
  }
  else
    mooseError("InertialForce: Either all or none of `beta`, `gamma`, `velocity`and `acceleration` "
               "should be provided as input.");

  // Check if HHT and explicit are being used simultaneously
  if (_alpha != 0 && _time_integrator.isExplicit())
    mooseError("InertialForce: HHT time integration parameter can only be used with Newmark-Beta "
               "time integration.");

  // Check if beta and explicit are being used simultaneously
  if (_has_beta && _time_integrator.isExplicit())
    mooseError("InertialForce: Newmark-beta integration parameter, beta, cannot be provided along "
               "with an explicit time "
               "integrator.");

  if (!(_time_integrator.isLumped() && _time_integrator.isExplicit()))
    if (parameters.isParamSetByUser("density_scaling"))
      this->paramError(
          "density_scaling",
          "Density (mass) scaling can only be used in lumped mass, explicit simulations");

    
}

Real
InertialForceFluid::computeQpResidual()
{

  Real dens = 0.0;
  for (unsigned ph = 0; ph < _num_phases; ++ph)
  {
     dens +=  _fluid_density_qp[_qp][ph];
  }

  if (_dt == 0)
    return 0;
  else if (_has_beta)
  {
    auto accel = 1.0 / _beta *
                 (((_u[_qp] - (*_u_old)[_qp]) / (_dt * _dt)) - (*_vel_old)[_qp] / _dt -
                  (*_accel_old)[_qp] * (0.5 - _beta));
    auto vel =
        (*_vel_old)[_qp] + (_dt * (1.0 - _gamma)) * (*_accel_old)[_qp] + _gamma * _dt * accel;
    return _test[_i][_qp]  * dens *
           (accel + vel * _eta * (1.0 + _alpha) - _alpha * _eta * (*_vel_old)[_qp]);
  }
  // Consistent mass option
  // Same for explicit, implicit, and implicit with HHT
  else
    return _test[_i][_qp]  * dens /_porosity[_qp]*
           ((*_u_dotdot_factor)[_qp] + (*_u_dot_factor)[_qp] * _eta * (1.0 + _alpha) -
            _alpha * _eta * (*_u_dot_old)[_qp]);
}


Real
InertialForceFluid::computeQpJacobian()
{
  Real dens = 0.0;
  for (unsigned ph = 0; ph < _num_phases; ++ph)
  {
     dens +=  _fluid_density_qp[_qp][ph];
  }

  if (_dt == 0)
    return 0;
  else

  {
    if (_has_beta)
      return _test[_i][_qp] * dens /_porosity[_qp] / (_beta * _dt * _dt) * _phi[this->_j][_qp] +
             _eta * (1 + _alpha) * _test[_i][_qp] * dens * _gamma / _beta / _dt *
                 _phi[this->_j][_qp];
    else
      return _test[_i][_qp] * dens /_porosity[_qp] * (*_du_dotdot_du)[_qp] * _phi[this->_j][_qp] +
             _eta * (1 + _alpha) * _test[_i][_qp] * dens  * (*_du_dot_du)[_qp] *
                 _phi[this->_j][_qp];
  }
}
