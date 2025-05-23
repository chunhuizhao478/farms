#include "AnalyticalSolution.h"
#include <cmath>

registerMooseObject("farmsApp", AnalyticalSolution);

InputParameters
AnalyticalSolution::validParams()
{
  InputParameters params = Function::validParams();
  params.addClassDescription("Analytical solution for MMS verification");
  params.addParam<Real>("amplitude", 1.0e-3, "Amplitude of the solution (m)");
  params.addParam<Real>("sigma", 50.0, "Spatial width parameter (m)");
  params.addParam<Real>("x0", 0.0, "x-coordinate of center (m)");
  params.addParam<Real>("y0", 0.0, "y-coordinate of center (m)");
  params.addParam<Real>("t0", 0.05, "Time of peak amplitude (s)");
  params.addParam<Real>("t_width", 0.02, "Temporal width of pulse (s)");
  params.addParam<Real>("omega", 2.0*M_PI*5.0, "Angular frequency (rad/s)");
  params.addParam<unsigned int>("component", 0, "Component: 0=x, 1=y");
  return params;
}

AnalyticalSolution::AnalyticalSolution(const InputParameters & parameters) :
  Function(parameters),
  _amplitude(getParam<Real>("amplitude")),
  _sigma(getParam<Real>("sigma")),
  _x0(getParam<Real>("x0")),
  _y0(getParam<Real>("y0")),
  _t0(getParam<Real>("t0")),
  _t_width(getParam<Real>("t_width")),
  _omega(getParam<Real>("omega")),
  _component(getParam<unsigned int>("component"))
{
}

Real
AnalyticalSolution::value(Real t, const Point & p) const
{
  // Get coordinates
  Real x = p(0);
  Real y = p(1);
  
  // Compute analytical solution
  Real r2 = (x - _x0)*(x - _x0) + (y - _y0)*(y - _y0);
  
  // Base function for both components
  Real u_base = _amplitude * std::exp(-r2/(2.0*_sigma*_sigma)) * 
                std::exp(-(t - _t0)*(t - _t0)/(2.0*_t_width*_t_width)) * 
                std::sin(_omega * t);
  
  // Direction factor based on component
  if (_component == 0) { // x-component
    if (std::abs(x - _x0) < 1e-10) {
      return 0.0; // To avoid division by zero
    }
    return u_base * (x - _x0) / std::sqrt(r2);
  }
  else { // y-component
    if (std::abs(y - _y0) < 1e-10) {
      return 0.0; // To avoid division by zero
    }
    return u_base * (y - _y0) / std::sqrt(r2);
  }
}