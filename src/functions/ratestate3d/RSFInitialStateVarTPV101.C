/*
Define Function for Spatially Variable Initial State Variable
Problem-Specific: TPV101-3D

According to SCEC TPV101 benchmark (Eq. 6):
θ_ini(x,y) = (L/V_0) * exp[(a*ln(2*sinh(τ_ini/(a*σ_ini))) - f_0 - a(x,y)*ln(V_ini/V_0)) / b]

The initial state variable must vary spatially to maintain uniform initial
shear stress given the spatially variable 'a' parameter.
*/

#include "RSFInitialStateVarTPV101.h"
#include <cmath>

registerMooseObject("farmsApp", RSFInitialStateVarTPV101);

InputParameters
RSFInitialStateVarTPV101::validParams()
{
  InputParameters params = Function::validParams();
  params.addClassDescription("Spatially variable initial state variable for TPV101 benchmark.");

  // RSF parameters
  params.addParam<Real>("f_0", 0.6, "Reference friction coefficient");
  params.addParam<Real>("V_0", 1e-6, "Reference slip velocity (m/s)");
  params.addParam<Real>("a_0", 0.008, "Base value of a in velocity-weakening region");
  params.addParam<Real>("b", 0.012, "Evolution effect parameter");
  params.addParam<Real>("L", 0.02, "Characteristic slip distance (m)");
  params.addParam<Real>("delta_a_0", 0.008, "Maximum increase in a for velocity-strengthening");

  // Initial conditions
  params.addParam<Real>("tau_ini", 75e6, "Initial shear stress (Pa)");
  params.addParam<Real>("sigma_ini", 120e6, "Initial normal stress (Pa)");
  params.addParam<Real>("V_ini", 1e-12, "Initial slip velocity (m/s)");

  // Geometry parameters
  params.addParam<Real>("W", 15000.0, "Half-width of velocity-weakening region (m)");
  params.addParam<Real>("w", 3000.0, "Transition layer width (m)");
  params.addParam<Real>("y_0", 7500.0, "Depth of center of VW region (m), positive value");

  return params;
}

RSFInitialStateVarTPV101::RSFInitialStateVarTPV101(const InputParameters & parameters)
  : Function(parameters),
    _f_0(getParam<Real>("f_0")),
    _V_0(getParam<Real>("V_0")),
    _a_0(getParam<Real>("a_0")),
    _b(getParam<Real>("b")),
    _L(getParam<Real>("L")),
    _delta_a_0(getParam<Real>("delta_a_0")),
    _tau_ini(getParam<Real>("tau_ini")),
    _sigma_ini(getParam<Real>("sigma_ini")),
    _V_ini(getParam<Real>("V_ini")),
    _W(getParam<Real>("W")),
    _w(getParam<Real>("w")),
    _y_0(getParam<Real>("y_0"))
{
}

Real
RSFInitialStateVarTPV101::boxcarB(Real x, Real W, Real w) const
{
  Real abs_x = std::abs(x);

  if (abs_x <= W)
  {
    return 1.0;
  }
  else if (abs_x >= W + w)
  {
    return 0.0;
  }
  else
  {
    // Transition region: W < |x| < W + w
    Real term1 = w / (abs_x - W - w);
    Real term2 = w / (abs_x - W);
    return 0.5 * (1.0 + std::tanh(term1 + term2));
  }
}

Real
RSFInitialStateVarTPV101::computeLocalA(const Point & p) const
{
  Real x_coord = p(0);
  Real z_coord = p(2);
  Real depth = std::abs(z_coord);

  Real B_x = boxcarB(x_coord, _W, _w);
  Real B_y = boxcarB(depth - _y_0, _W / 2.0, _w);

  Real delta_a = _delta_a_0 * (1.0 - B_x * B_y);

  return _a_0 + delta_a;
}

Real
RSFInitialStateVarTPV101::value(Real /*t*/, const Point & p) const
{
  // Get local a(x,y) value
  Real a_local = computeLocalA(p);

  // Compute initial state variable according to Eq. 6:
  // θ_ini = (L/V_0) * exp[(a*ln(2*sinh(τ_ini/(a*σ_ini))) - f_0 - a*ln(V_ini/V_0)) / b]

  // Step 1: Compute τ_ini / (a * σ_ini)
  Real tau_over_a_sigma = _tau_ini / (a_local * _sigma_ini);

  // Step 2: Compute 2*sinh(τ_ini/(a*σ_ini))
  // For large arguments, sinh(x) ≈ exp(x)/2, so 2*sinh(x) ≈ exp(x)
  // But we compute it directly for accuracy
  Real two_sinh = 2.0 * std::sinh(tau_over_a_sigma);

  // Step 3: Compute a * ln(2*sinh(...))
  Real term1 = a_local * std::log(two_sinh);

  // Step 4: Compute a * ln(V_ini / V_0)
  Real term2 = a_local * std::log(_V_ini / _V_0);

  // Step 5: Compute the exponent: [term1 - f_0 - term2] / b
  Real exponent = (term1 - _f_0 - term2) / _b;

  // Step 6: Compute θ_ini = (L/V_0) * exp(exponent)
  Real theta_ini = (_L / _V_0) * std::exp(exponent);

  return theta_ini;
}
