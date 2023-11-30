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

  //Parameters
  //Define pi
  Real pi = 3.14159265358979323846;
  //Define tolerance
  Real epsilon = 1e-12;
  //Define num_intervals
  int numofintervals = 1000;

  //compute R
  Real x_center = 0;
  Real y_center = 0;
  Real x_coord = p(0) - x_center; //along the strike direction
  Real y_coord = p(1) - y_center; //along the normal direction
  Real R = sqrt(x_coord*x_coord+y_coord*y_coord); //assume injection location is (0,0)

  //initialize pressure
  Real pressure = 0.0;

  //undrained lame constant
  Real drained_lambda   = 2 * _shear_modulus_mu *     _drained_nu / ( 1 - 2 *     _drained_nu );
  Real undrained_lambda = 2 * _shear_modulus_mu * _undrained_nu_u / ( 1 - 2 * _undrained_nu_u );

  //hydraulic diffusivity
  Real c = ( _permeability_k * ( undrained_lambda - drained_lambda ) * ( drained_lambda + 2 * _shear_modulus_mu ) ) / ( _viscosity_eta * _biotcoeff_alpha * _biotcoeff_alpha * ( undrained_lambda + 2 * _shear_modulus_mu ) );

  //Define z
  Real z = R * R / ( 4 * c * t );

  Real expIntz = 0.0;

  // Check if the argument is zero or negative
  if (z <= 0)
  {
    // Return an error message
    expIntz = 0.9 * 120e6;
  }
  else{

    // Set the upper limit of the integral to a large value
    double b = 10000;

    // Set the number of subintervals for the trapezoidal rule
    int n = 1000;

    // Compute the step size
    double h = (b - z) / n;

    // Initialize the sum
    double sum = 0;

    // Loop over the subintervals
    for (int i = 0; i < n; i++)
    {
      // Compute the endpoints of the subinterval
      double x1 = z + i * h;
      double x2 = z + (i + 1) * h;

      // Compute the function values at the endpoints
      double f1 = exp(-x1) / x1;
      double f2 = exp(-x2) / x2;

      // Add the area of the trapezoid to the sum
      sum += (f1 + f2) * h / 2;
    }

    expIntz = sum;

  }

  //compute pressure
  pressure = ( _flux_q * _viscosity_eta ) / ( 4 * pi * _density_rho_0 * _permeability_k ) * expIntz;

  return pressure;

}