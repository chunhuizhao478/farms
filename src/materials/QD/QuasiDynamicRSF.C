#include "QuasiDynamicRSF.h"
#include "InterfaceKernel.h"
#include "NestedSolve.h"
#include "libmesh/utility.h"
#include "RankTwoTensor.h"
#include "libmesh/vector_value.h"
#include "FEProblem.h"

registerMooseObject("farmsApp", QuasiDynamicRSF);

InputParameters
QuasiDynamicRSF::validParams()
{ 
  InputParameters params = CZMComputeLocalTractionTotalBaseQDRSF2D::validParams();
  params.addClassDescription("Quasi-Dynamic Rate-and-State Frictional Law.");
  
  // Rate-and-state friction parameters
  params.addRequiredParam<Real>("V_o", "Reference slip rate");
  params.addRequiredParam<Real>("f_o", "Initial friction coefficient");
  params.addRequiredParam<Real>("a", "Direct effect parameter");
  params.addRequiredParam<Real>("b", "State variable evolution parameter");
  params.addRequiredParam<Real>("L", "State variable characteristic distance");
  
  // Radiation damping parameters
  params.addRequiredParam<Real>("shear_modulus", "Shear modulus of the material");
  params.addRequiredParam<Real>("shear_wave_velocity", "Shear wave velocity");
  
  // Input variables
  params.addRequiredCoupledVar("normal_traction", "Normal traction input");
  params.addRequiredCoupledVar("tangential_traction", "Tangential traction input");
  params.addRequiredCoupledVar("pore_pressure", "Pore pressure input");
  
  // Optional parameters for enhanced weakening
  params.addParam<bool>("enhanced_weakening", false, "Flag for enhanced weakening");
  params.addParam<Real>("f_w", 0.0, "Weakened friction coefficient");
  params.addParam<Real>("V_w", 1.0, "Weakening slip rate");

  // Background stress parameters
  params.addParam<Real>("background_normal_stress", 0.0, "Background normal stress");
  params.addParam<Real>("background_tangential_stress", 0.0, "Background tangential stress");
  
  return params;
}

QuasiDynamicRSF::QuasiDynamicRSF(const InputParameters & parameters)
  : CZMComputeLocalTractionTotalBaseQDRSF2D(parameters),
    _V_o(getParam<Real>("V_o")),
    _f_o(getParam<Real>("f_o")),
    _a(getParam<Real>("a")),
    _b(getParam<Real>("b")),
    _L(getParam<Real>("L")),
    // Radiation damping
    _shear_modulus(getParam<Real>("shear_modulus")),
    _shear_wave_velocity(getParam<Real>("shear_wave_velocity")),
    _nu(_shear_modulus / (2.0 * _shear_wave_velocity)),
    // Enhanced weakening
    _enhanced_weakening(getParam<bool>("enhanced_weakening")),
    _f_w(getParam<Real>("f_w")),
    _V_w(getParam<Real>("V_w")),
    // Coupled input variables
    _normal_traction(coupledValue("normal_traction")),
    _tangential_traction(coupledValue("tangential_traction")),
    _p_p(coupledNeighborValue("pore_pressure")),
    _p_n(coupledValue("pore_pressure")),
    // Background stress
    _background_normal_stress(getParam<Real>("background_normal_stress")),
    _background_tangential_stress(getParam<Real>("background_tangential_stress")),
    // Old values
    _sliprate_old(getMaterialPropertyOldByName<Real>("sliprate")),
    _statevar_old(getMaterialPropertyOldByName<Real>("statevar")),
    _statevar_dot_old(getMaterialPropertyOldByName<Real>("statevar_dot")),
    _slip_old(getMaterialPropertyOldByName<Real>("slip"))
{
}

void
QuasiDynamicRSF::computeInterfaceTractionAndDerivatives()
{   
    // Get current inputs
    Real Tx = _background_tangential_stress + _tangential_traction[_qp];  // Tangential traction input
    Real Ty_input = _background_normal_stress - _normal_traction[_qp];  // Normal traction input
    Real p_f = std::max(_p_p[_qp],_p_n[_qp]);  // Pore pressure input
    
    // Compute effective normal traction
    Real Ty = Ty_input - p_f;

    // Get old state variables
    Real V_old = _sliprate_old[_qp];
    Real theta_old = _statevar_old[_qp];
    Real theta_dot_old = _statevar_dot_old[_qp];
    Real slip_old = _slip_old[_qp];

    // Slip rate solver
    Real V_new = solveSlipRate(Tx, Ty, p_f, theta_old, V_old);

    // Explicit time integration for slip
    Real slip_new = computeSlip(slip_old, V_new, V_old);

    // Update state variable (slip law)
    Real f_new = _a * std::asinh(V_new / (2.0 * _V_o) * std::exp(theta_old/ _a));
    Real theta_dot = - V_new / _L * (f_new - _f_o + (_b - _a) * std::log(V_new  / _V_o));
    Real theta_new = computeStateVariable(theta_old, theta_dot, theta_dot_old);

    // Store results
    _sliprate[_qp] = V_new;
    _statevar[_qp] = theta_new;
    _statevar_dot[_qp] = theta_dot;
    _slip[_qp] = slip_new;
}

Real 
QuasiDynamicRSF::solveSlipRate(Real Tx, Real Ty, Real p_avg, Real theta_old, Real V_old)
{
    // Initial guess
    Real x = V_old;
    
    // Check if slip rate is very small
    Real test = _V_o * std::exp(1.0 / _a * (Tx / (Ty) - theta_old));
    
    if (V_old < 1e-6)
    {
        return test;
    }

    // Newton-Raphson solver
    Real er = 1.0;
    Real x_n = 0.0;
    Real xo = V_old;
    Real xo_1 = xo - 0.0001 * xo;
    int iter = 1;
    int max_itr = 100000;

    while (er > 1e-8 && iter < max_itr)
    {
        Real Func, Func_1, dFunc;

        if (_enhanced_weakening)
        {
            // Enhanced weakening formulation
            Func = Tx - (_nu * xo) - 
                   (Ty) * (_f_w + 
                   (_a * std::asinh(xo / (2.0 * _V_o) * 
                    std::exp((_f_o - _f_w + _b * std::log(_V_o * theta_old / _L)) / _a))) / 
                   (1.0 + _L / (_V_w * theta_old)));

            Func_1 = Tx - (_nu * xo_1) - 
                     (Ty) * (_f_w + 
                     (_a * std::asinh(xo_1 / (2.0 * _V_o) * 
                      std::exp((_f_o - _f_w + _b * std::log(_V_o * theta_old / _L)) / _a))) / 
                     (1.0 + _L / (_V_w * theta_old)));
        }
        else
        {
            // Standard Dieterich-Ruina formulation
            Func = Tx - (_nu * xo) - 
                   (Ty) * (_a * std::asinh(xo / (2.0 * _V_o) * 
                    std::exp(theta_old/ _a)));

            Func_1 = Tx - (_nu * xo_1) - 
                     (Ty) * (_a * std::asinh(xo_1 / (2.0 * _V_o) * 
                      std::exp(theta_old/ _a)));
        }

        // Update slip rate
        if (er > 1.0)
        {
            // Secant method for first iterations
            x_n = xo - Func * (xo - xo_1) / (Func - Func_1);
            er = std::abs(x_n - xo) / std::abs(x_n);
        }
        else
        {
            // Newton-Raphson method
            if (_enhanced_weakening)
            {
                dFunc = -_nu / 2.0 - 
                        (Ty * _a * std::exp((_f_o - _f_w + _b * std::log(_V_o * theta_old / _L)) / _a)) / 
                        (2.0 * _V_o * (1.0 + _L / (_V_w * theta_old)) * 
                         std::sqrt(1.0 + std::pow(xo * std::exp((_f_o - _f_w + _b * std::log(_V_o * theta_old / _L)) / _a) / (2.0 * _V_o), 2)));
            }
            else
            {
                dFunc = -_nu / 2.0 - 
                        (Ty * _a * std::exp(theta_old/ _a)) / 
                        (2.0 * _V_o * 
                         std::sqrt(1.0 + std::pow(xo * std::exp(theta_old/ _a) / (2.0 * _V_o), 2)));
            }

            x_n = xo - Func / dFunc;
            er = std::abs(x_n - xo) / std::abs(x_n);
        }

        // Update for next iteration
        xo_1 = xo;
        xo = x_n;
        iter++;
    }

    if (iter >= max_itr)
    {
        mooseWarning("Quasi-Dynamic RSF: Slip rate solver did not converge");
    }

    return std::abs(x_n);
}

Real 
QuasiDynamicRSF::computeSlip(Real slip_old, Real V_new, Real V_old)
{
    // Explicit time integration for slip using Adams-Bashforth-like method
    return slip_old + _dt/2.0 * (3.0 * V_new - V_old);
}

Real 
QuasiDynamicRSF::computeStateVariable(Real theta_old, Real theta_dot, Real theta_dot_old)
{
    // Explicit time integration for state variable
    return theta_old + _dt/2.0 * (3.0 * theta_dot - theta_dot_old);
}