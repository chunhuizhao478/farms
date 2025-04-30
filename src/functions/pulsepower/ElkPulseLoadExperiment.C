#include "ElkPulseLoadExperiment.h"

registerMooseObject("farmsApp", ElkPulseLoadExperiment);

InputParameters
ElkPulseLoadExperiment::validParams()
{
  InputParameters params = Function::validParams();
  params.addRequiredParam<Real>("shape_param_alpha","shape parameter alpha");
  params.addRequiredParam<Real>("shape_param_beta","shape parameter beta");
  params.addRequiredParam<Real>("rise_time","rise time t0 (second)");
  params.addRequiredParam<Real>("single_pulse_duration","single pulse duration (s)");
  params.addRequiredParam<Real>("EM","energy of the pulse (kJ)");
  params.addRequiredParam<Real>("gap","the gap distance of which the discharge takes place");
  params.addRequiredParam<Real>("convert_efficiency","convert efficiency eta (EB = eta EM)");
  params.addRequiredParam<Real>("fitting_param_alpha","fitting parameter for peak pressure");
  params.addRequiredParam<std::vector<Real>>("discharge_center", "discharge center (x,y,z) (m)");
  params.addRequiredParam<int>("number_of_pulses","number of pulse, assume pulses are continuous");
  params.addParam<Real>("peak_pressure", 0.0, "The peak value of pressure, if not specified, it will be calculated based on the Dsensor.");
  params.addParam<Real>("minimum_applied_pressure", 0.0, "Minimum applied pressure to mimic the effect of water pressure");
  params.addParam<bool>("use_minimum_applied_pressure", false, "Flag to use minimum applied pressure");
  return params;
}

ElkPulseLoadExperiment::ElkPulseLoadExperiment(const InputParameters & parameters)
  : Function(parameters),
  _shape_param_alpha(getParam<Real>("shape_param_alpha")),
  _shape_param_beta(getParam<Real>("shape_param_beta")),
  _rise_time(getParam<Real>("rise_time")),
  _single_pulse_duration(getParam<Real>("single_pulse_duration")),
  _convert_efficiency(getParam<Real>("convert_efficiency")),
  _EM(getParam<Real>("EM")),
  _gap(getParam<Real>("gap")),
  _fitting_param_alpha(getParam<Real>("fitting_param_alpha")),
  _discharge_center(getParam<std::vector<Real>>("discharge_center")),
  _number_of_pulses(getParam<int>("number_of_pulses")),
  _peak_pressure(getParam<Real>("peak_pressure")),
  _minimum_applied_pressure(getParam<Real>("minimum_applied_pressure")),
  _use_minimum_applied_pressure(getParam<bool>("use_minimum_applied_pressure")) 
{
  if (_use_minimum_applied_pressure == true && _use_minimum_applied_pressure < 0.0){
    mooseError("Minimum applied pressure must be non-negative");
  }
}

Real
ElkPulseLoadExperiment::value(Real t, const Point & p) const
{

  // Peak stress Pp
  
  // Get coordinate
  Real xcoord = p(0); //along the x direction
  Real ycoord = p(1); //along the y direction
  Real zcoord = p(2); //along the z direction

  Real Dsensor = 0.0;
  // Within the gap, assume uniform Dsensor
  if ( (zcoord >= _discharge_center[2] - _gap/2.0) and (zcoord <= _discharge_center[2] + _gap/2.0)){
    Dsensor = (1000) * std::sqrt(std::pow(xcoord-_discharge_center[0],2)+std::pow(ycoord-_discharge_center[1],2));
  }
  else if( zcoord < _discharge_center[2] - _gap/2.0 ){
    Dsensor = (1000) * std::sqrt(std::pow(xcoord-_discharge_center[0],2)+std::pow(ycoord-_discharge_center[1],2)+std::pow(zcoord-(_discharge_center[2] - _gap/2.0),2));
  }
  else if( zcoord > _discharge_center[2] + _gap/2.0 ){
    Dsensor = (1000) * std::sqrt(std::pow(xcoord-_discharge_center[0],2)+std::pow(ycoord-_discharge_center[1],2)+std::pow(zcoord-(_discharge_center[2] + _gap/2.0),2));
  }

  // Estimate peak pressure (bar mm KJ) -> Pp (bar -> Pa)
  Real Pp = (0.1 * 1e6) * 9000 * 1.0 / Dsensor * std::pow(_convert_efficiency*_EM,_fitting_param_alpha);

  // Check if the peak pressure is specified
  if (_peak_pressure > 0.0){
    Pp = _peak_pressure;
  }

  // Constants
  Real PEAK_MAGNITUDE = Pp;                          // Peak magnitude in Pascals
  Real PULSE_DURATION_US = _single_pulse_duration;   // Duration of one pulse in microseconds
  int TOTAL_PULSES = _number_of_pulses;              // Total number of pulses

  //calculate the magnitude at a given time
  Real mod_time = std::fmod(t, PULSE_DURATION_US); // Time within the current pulse period

  Real total_duration_s = TOTAL_PULSES * PULSE_DURATION_US; // Total duration in seconds
  
  Real pulse_load = 0.0;
  
  if ( t <= total_duration_s ){
    pulse_load = PEAK_MAGNITUDE * (std::exp(-_shape_param_alpha * mod_time) - std::exp(-_shape_param_beta * mod_time)) / 
                  (std::exp(-_shape_param_alpha * _rise_time) - std::exp(-_shape_param_beta * _rise_time));
  }
  else{
    pulse_load = 0.0;
  }

  // Apply minimum applied pressure
  if (_use_minimum_applied_pressure){
    pulse_load = std::max(pulse_load, _minimum_applied_pressure);
  }

  return pulse_load;

}