#include "MMSSourceTermY.h"
#include <cmath>

registerMooseObject("farmsApp", MMSSourceTermY);

InputParameters
MMSSourceTermY::validParams()
{
  InputParameters params = BodyForce::validParams();
  params.addClassDescription("Implements the source term for the MMS in y-direction (continuous domain)");
  params.addRequiredParam<Real>("delta", "Total displacement amplitude (m)");
  params.addRequiredParam<Real>("R_m", "Characteristic length (m)");
  params.addRequiredParam<Real>("tbar", "Rupture time (s)");
  params.addRequiredParam<Real>("tw", "Event duration parameter (s)");
  params.addRequiredParam<Real>("Vmin", "Minimum slip velocity (m/s)");
  params.addRequiredParam<Real>("density", "the density");
  params.addRequiredParam<Real>("lambda", "lame's constant");
  params.addRequiredParam<Real>("shear_modulus", "Shear modulus");
  params.addParam<Real>("tau_background", 40e6, "Background stress (Pa)");
  return params;
}

MMSSourceTermY::MMSSourceTermY(const InputParameters & parameters) :
  BodyForce(parameters),
  _delta(getParam<Real>("delta")),
  _R_m(getParam<Real>("R_m")),
  _tbar(getParam<Real>("tbar")),
  _tw(getParam<Real>("tw")),
  _Vmin(getParam<Real>("Vmin")),
  _density(getParam<Real>("density")),
  _lambda(getParam<Real>("lambda")),
  _mu(getParam<Real>("shear_modulus")),
  _tau_background(getParam<Real>("tau_background"))
{
}

Real
MMSSourceTermY::computeQpResidual()
{
  // Get the coordinates and time
  Real x = _q_point[_qp](0);
  Real y = _q_point[_qp](1);
  Real t = _t + _dt;
  
  // Calculate the source term
  Real source = sourceTerm(x, y, t);
  
  // Return negative of source term (MOOSE convention)
  return -source * _test[_i][_qp];
}

Real
MMSSourceTermY::spatialFunction(Real x, Real y) const
{
  // Continuous exponential function without sign change
  Real width = _R_m / 2.0;
  return std::exp(-(x*x + y*y)/(2.0*width*width));
}

Real
MMSSourceTermY::temporalFunction(Real t) const
{
    return (_delta/2.0) * std::exp(( t - _tbar)/_tw);
}

Real
MMSSourceTermY::temporalDerivative(Real t) const
{
    return (_delta/2.0) * std::exp(( t - _tbar)/_tw)/_tw;
}

Real
MMSSourceTermY::temporalSecondDerivative(Real t) const
{
  
    return (_delta/2.0) * std::exp(( t - _tbar)/_tw)/_tw/_tw;
}

Real
MMSSourceTermY::sourceTerm(Real x, Real y, Real t) const
{
Real width = _R_m / 2.0;
  Real spatial_phi = std::exp(-(x*x + y*y)/(2.0*width*width));
  Real d_spatial_phi_dx = -x/(width*width) * spatial_phi;
  Real d_spatial_phi_dy = -y/(width*width) * spatial_phi;
  
  // Second spatial derivatives
  Real d2_spatial_phi_dx2 = ((x*x)/(width*width*width*width) - 1.0/(width*width)) * spatial_phi;
  Real d2_spatial_phi_dy2 = ((y*y)/(width*width*width*width) - 1.0/(width*width)) * spatial_phi;
  Real d2x_spatial_phi_dy2 = ((x*y)/(width*width*width*width)) * spatial_phi;
  Real d2y_spatial_phi_dx2 = ((x*y)/(width*width*width*width)) * spatial_phi;
  
  // Temporal functions
  Real K = temporalFunction(t);
  Real dK_dt = temporalDerivative(t);
  Real d2K_dt2 = temporalSecondDerivative(t);
  
  // Displacement components
  Real displacement = K * spatial_phi;
  Real disp_dx = K * d_spatial_phi_dx;
  Real disp_dy = K * d_spatial_phi_dy;
  
  // Strain components (small strain assumption)
  Real epsilon_xx = disp_dx;
  Real epsilon_yy = disp_dy;
  Real epsilon_xy = 0.5 * (disp_dy + disp_dx);
  
  // Stress components (Hooke's law)
  Real sigma_xx = (_lambda + 2.0*_mu) * epsilon_xx;
  Real sigma_yy = (_lambda + 2.0*_mu) * epsilon_yy;
  Real sigma_xy = _mu * 2.0 * epsilon_xy;
  
  // Stress divergence
  Real div_sigma_y = 
    _lambda * K * (d2_spatial_phi_dx2 + d2_spatial_phi_dy2) +  2.0 * _mu * K * d2_spatial_phi_dy2;
   
  Real div_sigma_xy = 
       _mu * K * 1/2* (d2x_spatial_phi_dy2+d2y_spatial_phi_dx2);
  
  // Acceleration term
  Real accel_term = _density * d2K_dt2 * spatial_phi;
  
  // Final source term
  return accel_term - div_sigma_y - div_sigma_xy;
}