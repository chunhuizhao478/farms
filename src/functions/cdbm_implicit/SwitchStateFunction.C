#include "SwitchStateFunction.h"

registerMooseObject("farmsApp", SwitchStateFunction);

InputParameters
SwitchStateFunction::validParams()
{
  InputParameters params = Function::validParams();
  params.addRequiredParam<PostprocessorName>("vel_sol", "Postprocessor value for velocity solution");
  params.addRequiredParam<Real>("vel_threshold_quasi_to_dyn", "Velocity threshold from quasi-dynamic to dynamic");
  params.addRequiredParam<Real>("vel_threshold_dyn_to_quasi", "Velocity threshold from dynamic to quasi-dynamic");
  params.addRequiredParam<PostprocessorName>("stateflag_old", "Flag for identifying the state of the simulation (Dynamic or Quasi-dynamic)");
  return params;
}

SwitchStateFunction::SwitchStateFunction(const InputParameters & parameters)
  : Function(parameters),
    _vel_sol(getPostprocessorValue("vel_sol")),
    _vel_threshold_quasi_to_dyn(getParam<Real>("vel_threshold_quasi_to_dyn")),
    _vel_threshold_dyn_to_quasi(getParam<Real>("vel_threshold_dyn_to_quasi")),
    _stateflag_old(getPostprocessorValue("stateflag_old")),
    _stateflag_older(getPostprocessorValueOld("stateflag_old")),
    _stateflag_oldolder(getPostprocessorValueOlder("stateflag_old"))
{
}

Real
SwitchStateFunction::value(Real t, const Point & /*p*/) const
{

  Real state = 0;

  //skip the initial inactive step //we actually don't need this, but keeping it for now
  //the condition will be skipped in the first time step, see FarmsConditionalEnableControl
  if (t < _dt){
    return 0; //do nothing
  }

  //if the velocity solution is less than the threshold, current state is dynamic, switch to quasi-dynamic
  //the state has been dynamic for the last 3 time steps
  if ((_vel_sol < _vel_threshold_dyn_to_quasi) && (_stateflag_old == 1) && (_stateflag_older == 1) && (_stateflag_oldolder == 1)){
    state = 0; //switch to quasi-dynamic state
  }
  //if the velocity solution is greater than the threshold, current state is quasi-dynamic, switch to dynamic
  //the state has been quasi-dynamic for the last 3 time steps
  else if ((_vel_sol > _vel_threshold_quasi_to_dyn) && (_stateflag_old == 0) && (_stateflag_older == 0) && (_stateflag_oldolder == 0)){
    state = 1; //switch to dynamic state
  }
  //if the velocity solution is between the thresholds, keep the current state
  else{
    state = _stateflag_old;
  }

  return state;

}