#include "GaussianSourceTerm.h"
#include <cmath>

registerMooseObject("farmsApp", GaussianSourceTerm);

InputParameters
GaussianSourceTerm::validParams()
{
  InputParameters params = BodyForce::validParams();
  params.addClassDescription("Simple Gaussian pulse source for elastic wave propagation with MMS");
  params.addParam<Real>("amplitude", 1.0e-3, "Amplitude of the source (m)");
  params.addParam<Real>("sigma", 50.0, "Spatial width parameter (m)");
  params.addParam<Real>("x0", 0.0, "x-coordinate of center (m)");
  params.addParam<Real>("y0", 0.0, "y-coordinate of center (m)");
  params.addParam<Real>("t0", 0.05, "Time of peak amplitude (s)");
  params.addParam<Real>("t_width", 0.02, "Temporal width of pulse (s)");
  params.addParam<Real>("omega", 2.0*M_PI*5.0, "Angular frequency (rad/s)");
  params.addParam<unsigned int>("component", 0, "Component: 0=x, 1=y");
  return params;
}

GaussianSourceTerm::GaussianSourceTerm(const InputParameters & parameters) :
  BodyForce(parameters),
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
GaussianSourceTerm::computeQpResidual()
{
  // Get coordinates and time
  Real x = _q_point[_qp](0);
  Real y = _q_point[_qp](1);
  Real t = _t;
  
  // Calculate MMS source term
  // For MMS, we need to compute the divergence of stress and acceleration terms
  
  // First, calculate the analytical solution
  Real r2 = (x - _x0)*(x - _x0) + (y - _y0)*(y - _y0);
  Real u_mms = _amplitude * std::exp(-r2/(2.0*_sigma*_sigma)) * 
               std::exp(-(t - _t0)*(t - _t0)/(2.0*_t_width*_t_width)) * 
               std::sin(_omega * t);
  
  // For the MMS source term, we need:
  // f = rho * d²u/dt² - div(sigma)
  
  // Acceleration term: rho * d²u/dt²
  Real density = 8000.0; // kg/m³
  Real d2u_dt2 = _amplitude * std::exp(-r2/(2.0*_sigma*_sigma)) * 
                (
                  // Second time derivative of Gaussian
                  std::exp(-(t - _t0)*(t - _t0)/(2.0*_t_width*_t_width)) * 
                  (
                    ((t - _t0)*(t - _t0)/(_t_width*_t_width*_t_width*_t_width) - 1.0/(_t_width*_t_width))
                  ) * std::sin(_omega * t) +
                  
                  // First time derivative of Gaussian * first time derivative of sine
                  2.0 * (-(t - _t0)/(_t_width*_t_width)) * 
                  std::exp(-(t - _t0)*(t - _t0)/(2.0*_t_width*_t_width)) * 
                  _omega * std::cos(_omega * t) +
                  
                  // Gaussian * second time derivative of sine
                  std::exp(-(t - _t0)*(t - _t0)/(2.0*_t_width*_t_width)) * 
                  (-_omega*_omega) * std::sin(_omega * t)
                );
  
  Real accel_term = density * d2u_dt2;
  
  // Divergence of stress term: div(sigma)
  // For isotropic linear elasticity with Lame parameters lambda and mu
  Real lambda = 1.15e11; // Pa (for steel)
  Real mu = 7.69e10;     // Pa (for steel)
  
  // Compute spatial derivatives
  Real d2u_dx2 = _amplitude * 
                 (
                   // Second x-derivative of Gaussian
                   ((x - _x0)*(x - _x0)/(_sigma*_sigma*_sigma*_sigma) - 1.0/(_sigma*_sigma)) * 
                   std::exp(-r2/(2.0*_sigma*_sigma))
                 ) * 
                 std::exp(-(t - _t0)*(t - _t0)/(2.0*_t_width*_t_width)) * 
                 std::sin(_omega * t);
  
  Real d2u_dy2 = _amplitude * 
                 (
                   // Second y-derivative of Gaussian
                   ((y - _y0)*(y - _y0)/(_sigma*_sigma*_sigma*_sigma) - 1.0/(_sigma*_sigma)) * 
                   std::exp(-r2/(2.0*_sigma*_sigma))
                 ) * 
                 std::exp(-(t - _t0)*(t - _t0)/(2.0*_t_width*_t_width)) * 
                 std::sin(_omega * t);
  
  // Divergence depends on component
  Real div_stress = 0.0;
  
  if (_component == 0) { // x-component
    div_stress = (lambda + 2.0*mu) * d2u_dx2 + mu * d2u_dy2;
  }
  else { // y-component
    div_stress = (lambda + 2.0*mu) * d2u_dy2 + mu * d2u_dx2;
  }
  
  // Final source term
  Real source = accel_term - div_stress;
  
  // Return negative (MOOSE convention)
  return -source * _test[_i][_qp];
}