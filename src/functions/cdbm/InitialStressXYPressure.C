/* Reference: John W. RUDNICKI : FLUID MASS SOURCES AND POINT FORCES IN LINEAR ELASTIC DIFFUSIVE SOLIDS */

#include "InitialStressXYPressure.h"

registerMooseObject("farmsApp", InitialStressXYPressure);

InputParameters
InitialStressXYPressure::validParams()
{
  InputParameters params = Function::validParams();
  params.addRequiredParam<Real>(          "flux_q", "flux");
  params.addRequiredParam<Real>(   "density_rho_0", "fluid density");
  params.addRequiredParam<Real>(  "permeability_k", "permeability");
  params.addRequiredParam<Real>(   "viscosity_eta", "viscosity");
  params.addRequiredParam<Real>( "biotcoeff_alpha", "biot coefficient");
  params.addRequiredParam<Real>(  "undrained_nu_u", "undrained poisson's ratio");
  params.addRequiredParam<Real>("shear_modulus_mu", "shear modulus");
  params.addRequiredParam<Real>(      "drained_nu", "drained poisson's ratio");
  return params;
}

InitialStressXYPressure::InitialStressXYPressure(const InputParameters & parameters)
  : Function(parameters),
  _flux_q(getParam<Real>("flux_q")),
  _density_rho_0(getParam<Real>("density_rho_0")),
  _permeability_k(getParam<Real>("permeability_k")),
  _viscosity_eta(getParam<Real>("viscosity_eta")),
  _biotcoeff_alpha(getParam<Real>("biotcoeff_alpha")),
  _undrained_nu_u(getParam<Real>("undrained_nu_u")),
  _shear_modulus_mu(getParam<Real>("shear_modulus_mu")),
  _drained_nu(getParam<Real>("drained_nu"))
{
}

Real
InitialStressXYPressure::value(Real t, const Point & p) const
{

  //compute R
  Real x_coord = p(0); //along the strike direction
  Real y_coord = p(1); //along the normal direction
  Real R = sqrt(x_coord*x_coord+y_coord*y_coord); //assume injection location is (0,0)

  //initialize pressure
  Real pressure = 0.0;

  //undrained lame constant
  Real drained_lambda   = 2 * _shear_modulus_mu *     _drained_nu / ( 1 - 2 *     _drained_nu );
  Real undrained_lambda = 2 * _shear_modulus_mu * _undrained_nu_u / ( 1 - 2 * _undrained_nu_u );

  //hydraulic diffusivity
  Real c = ( _permeability_k * ( undrained_lambda - drained_lambda ) * ( drained_lambda + 2 * _shear_modulus_mu ) ) / ( _viscosity_eta * _biotcoeff_alpha * _biotcoeff_alpha * ( undrained_lambda + 2 * _shear_modulus_mu ) );

  //Euler's constant
  Real gamma = 0.57721;

  //Define pi
  Real pi = 3.14159265358979323846;

  //compute pressure
  pressure = ( _flux_q * _viscosity_eta ) / ( 4 * pi * _density_rho_0 * _permeability_k ) * ( -1 * gamma - log( R * R / ( 4 * c * t ) ) );

  return pressure;

}