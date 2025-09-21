#include "ElkPulseLoadExperimentWu2022Paper.h"

registerMooseObject("farmsApp", ElkPulseLoadExperimentWu2022Paper);

InputParameters
ElkPulseLoadExperimentWu2022Paper::validParams()
{
  InputParameters params = Function::validParams();
  params.addRequiredParam<Real>("shape_param_alpha","shape parameter alpha");
  params.addRequiredParam<Real>("shape_param_beta","shape parameter beta");
  params.addRequiredParam<Real>("rise_time","rise time t0 (second)");
  params.addRequiredParam<Real>("single_pulse_duration","single pulse duration (s)");
  params.addRequiredParam<std::vector<Real>>("discharge_center", "discharge center (x,y,z) (m)");
  params.addRequiredParam<std::vector<Real>>("Pmax_coefficients", "coefficients for Pmax");
  params.addRequiredParam<int>("number_of_pulses","number of pulse, assume pulses are continuous");
  params.addParam<Real>("minimum_applied_pressure", 0.0, "Minimum applied pressure to mimic the effect of water pressure");
  params.addParam<bool>("use_minimum_applied_pressure", false, "Flag to use minimum applied pressure");
  return params;
}

ElkPulseLoadExperimentWu2022Paper::ElkPulseLoadExperimentWu2022Paper(const InputParameters & parameters)
  : Function(parameters),
  _shape_param_alpha(getParam<Real>("shape_param_alpha")),
  _shape_param_beta(getParam<Real>("shape_param_beta")),
  _rise_time(getParam<Real>("rise_time")),
  _single_pulse_duration(getParam<Real>("single_pulse_duration")),
  _discharge_center(getParam<std::vector<Real>>("discharge_center")),
  _pmax_coefficients(getParam<std::vector<Real>>("Pmax_coefficients")),
  _number_of_pulses(getParam<int>("number_of_pulses")),
  _minimum_applied_pressure(getParam<Real>("minimum_applied_pressure")),
  _use_minimum_applied_pressure(getParam<bool>("use_minimum_applied_pressure")) 
{
  if (_use_minimum_applied_pressure == true && _use_minimum_applied_pressure < 0.0){
    mooseError("Minimum applied pressure must be non-negative");
  }
}

Real
ElkPulseLoadExperimentWu2022Paper::value(Real t, const Point & p) const
{

  // Peak stress Pp
  
  // Get coordinate
  Real xcoord = p(0); //along the x direction
  Real ycoord = p(1); //along the y direction
  Real zcoord = p(2); //along the z direction

  // Compute r relative to the discharge center
  Real r = std::sqrt(std::pow(xcoord - _discharge_center[0], 2) +
                     std::pow(ycoord - _discharge_center[1], 2) +
                     std::pow(zcoord - _discharge_center[2], 2));

  // Convert r from m to mm
  r *= 1000.0;

  // Get Pmax coefficients
  // Polynomial fit (degree 4) for r >= 1.6 mm:
  // Variables: r in mm, P in MPa
  // Coefficients: determined using pulseloading/pulse_load_singlepulse_wu2022.py
  // there are five coefficients
  Real a0 = _pmax_coefficients[0];
  Real a1 = _pmax_coefficients[1];
  Real a2 = _pmax_coefficients[2];
  Real a3 = _pmax_coefficients[3];
  Real a4 = _pmax_coefficients[4];

  Real Pmax = a0 * std::pow(r,4) + a1 * std::pow(r,3) + a2 * std::pow(r,2) + a3 * r + a4; //in MPa

  // Convert Pmax from MPa to Pa
  Real Pmax_Pa = Pmax * 1e6;

  // Constants
  Real PEAK_MAGNITUDE = Pmax_Pa;                     // Peak magnitude in Pascals
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