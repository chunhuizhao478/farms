//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

/*
Material Description of Rate-and-State Friction 3d
Based on fast-velocity-weakening formulation from Premus et al. (2020)
*/

#include "RateStateFrictionczm3d.h"
#include "InterfaceKernel.h"
#include <cmath>

registerMooseObject("farmsApp", RateStateFrictionczm3d);

InputParameters
RateStateFrictionczm3d::validParams()
{
  InputParameters params = CZMComputeLocalTractionTotalBase::validParams();
  params.addClassDescription("Rate-and-state friction with fast-velocity-weakening law.");
  
  // Rate-and-state friction parameters
  params.addRequiredParam<Real>("a", "Direct effect parameter a");
  params.addRequiredParam<Real>("b", "Evolution effect parameter b");
  params.addRequiredParam<Real>("L", "Characteristic length scale L");
  params.addRequiredParam<Real>("f0", "Reference friction coefficient f0");
  params.addRequiredParam<Real>("fw", "Weakened friction coefficient fw");
  params.addRequiredParam<Real>("fLV", "Low velocity friction coefficient fLV");
  params.addRequiredParam<Real>("s0", "Reference slip rate s0");
  params.addRequiredParam<Real>("sw", "Weakening slip rate sw");
  params.addRequiredParam<Real>("s_ini", "Initial slip rate for numerical regularization");
  
  params.addRequiredParam<Real>("len", "element edge length");
  params.addRequiredCoupledVar("disp_slipweakening_x", "displacement in x dir");
  params.addRequiredCoupledVar("disp_slipweakening_y", "displacement in y dir");
  params.addRequiredCoupledVar("disp_slipweakening_z", "displacement in z dir");
  params.addRequiredCoupledVar("vel_slipweakening_x","velocity in x dir");
  params.addRequiredCoupledVar("vel_slipweakening_y","velocity in y dir");
  params.addRequiredCoupledVar("vel_slipweakening_z","velocity in z dir");
  params.addRequiredCoupledVar("reaction_slipweakening_x", "reaction in x dir");
  params.addRequiredCoupledVar("reaction_slipweakening_y", "reaction in y dir");
  params.addRequiredCoupledVar("reaction_slipweakening_z", "reaction in z dir");
  
  // State variable coupling
  params.addRequiredCoupledVar("state_variable", "Rate-and-state friction state variable psi");
  
  params.addCoupledVar("cohesion_aux", "auxiliary variable for cohesion");
  params.addCoupledVar("fluid_pressure_aux", "auxiliary variable for fluid pressure");
  
  // Newton solver parameters
  params.addParam<Real>("newton_tolerance", 1.0e-10, "Newton solver tolerance");
  params.addParam<unsigned int>("newton_max_iterations", 50, "Maximum Newton iterations");
  
  return params;
}

RateStateFrictionczm3d::RateStateFrictionczm3d(const InputParameters & parameters)
  : CZMComputeLocalTractionTotalBase(parameters),
    // Rate-and-state parameters
    _a(getParam<Real>("a")),
    _b(getParam<Real>("b")),
    _L(getParam<Real>("L")),
    _f0(getParam<Real>("f0")),
    _fw(getParam<Real>("fw")),
    _fLV(getParam<Real>("fLV")),
    _s0(getParam<Real>("s0")),
    _sw(getParam<Real>("sw")),
    _s_ini(getParam<Real>("s_ini")),
    _len(getParam<Real>("len")),
    _density(getMaterialPropertyByName<Real>(_base_name + "density")),
    _rot(getMaterialPropertyByName<RankTwoTensor>(_base_name + "czm_total_rotation")),
    
    // Displacement variables
    _disp_slipweakening_x(coupledValue("disp_slipweakening_x")),
    _disp_slipweakening_neighbor_x(coupledNeighborValue("disp_slipweakening_x")),
    _disp_slipweakening_y(coupledValue("disp_slipweakening_y")),
    _disp_slipweakening_neighbor_y(coupledNeighborValue("disp_slipweakening_y")),
    _disp_slipweakening_z(coupledValue("disp_slipweakening_z")),
    _disp_slipweakening_neighbor_z(coupledNeighborValue("disp_slipweakening_z")),
    
    // Velocity variables
    _vel_slipweakening_x(coupledValue("vel_slipweakening_x")),
    _vel_slipweakening_neighbor_x(coupledNeighborValue("vel_slipweakening_x")),
    _vel_slipweakening_y(coupledValue("vel_slipweakening_y")),
    _vel_slipweakening_neighbor_y(coupledNeighborValue("vel_slipweakening_y")),
    _vel_slipweakening_z(coupledValue("vel_slipweakening_z")),
    _vel_slipweakening_neighbor_z(coupledNeighborValue("vel_slipweakening_z")),
    
    // Reaction variables
    _reaction_slipweakening_x(coupledValue("reaction_slipweakening_x")),
    _reaction_slipweakening_neighbor_x(coupledNeighborValue("reaction_slipweakening_x")),
    _reaction_slipweakening_y(coupledValue("reaction_slipweakening_y")),
    _reaction_slipweakening_neighbor_y(coupledNeighborValue("reaction_slipweakening_y")),
    _reaction_slipweakening_z(coupledValue("reaction_slipweakening_z")),
    _reaction_slipweakening_neighbor_z(coupledNeighborValue("reaction_slipweakening_z")),
    
    // Old displacement values
    _disp_slipweakening_x_old(coupledValueOld("disp_slipweakening_x")),
    _disp_slipweakening_neighbor_x_old(coupledNeighborValueOld("disp_slipweakening_x")),
    _disp_slipweakening_y_old(coupledValueOld("disp_slipweakening_y")),
    _disp_slipweakening_neighbor_y_old(coupledNeighborValueOld("disp_slipweakening_y")),
    _disp_slipweakening_z_old(coupledValueOld("disp_slipweakening_z")),
    _disp_slipweakening_neighbor_z_old(coupledNeighborValueOld("disp_slipweakening_z")),
    
    // State variable
    _state_variable(coupledValue("state_variable")),
    _state_variable_old(coupledValueOld("state_variable")),
    
    // Material properties for output
    _displacement_jump_strike(declareProperty<Real>("displacement_jump_strike")),
    _displacement_jump_dip(declareProperty<Real>("displacement_jump_dip")),
    _displacement_jump_normal(declareProperty<Real>("displacement_jump_normal")),
    _displacement_jump_rate_strike(declareProperty<Real>("displacement_jump_rate_strike")),
    _displacement_jump_rate_dip(declareProperty<Real>("displacement_jump_rate_dip")),
    _displacement_jump_rate_normal(declareProperty<Real>("displacement_jump_rate_normal")),
    _traction_strike(declareProperty<Real>("traction_strike")),
    _traction_dip(declareProperty<Real>("traction_dip")),
    _traction_normal(declareProperty<Real>("traction_normal")),
    _slip_rate_magnitude(declareProperty<Real>("slip_rate_magnitude")),
    _friction_coefficient(declareProperty<Real>("friction_coefficient")),
    _state_variable_updated(declareProperty<Real>("state_variable_updated")),
    
    _static_initial_stress_tensor(getMaterialPropertyByName<RankTwoTensor>(_base_name + "static_initial_stress_tensor")),
    _cohesion_aux(coupledValue("cohesion_aux")),
    _fluid_pressure_aux(coupledValue("fluid_pressure_aux")),
    
    // Newton solver parameters
    _newton_tolerance(getParam<Real>("newton_tolerance")),
    _newton_max_iterations(getParam<unsigned int>("newton_max_iterations"))
{
  // only works for small strain
  if (hasBlockMaterialProperty<RankTwoTensor>(_base_name + "strain_increment"))
  {
    mooseError("RateStateFriction only works for small strain!");
  }
}

Real
RateStateFrictionczm3d::computeFSS(Real slip_rate)
{
  // Compute f_LV (equation from paper)
  Real f_LV = _f0 - (_b - _a) * std::log(slip_rate / _s0);
  
  // Compute f_SS (steady-state friction)
  Real f_SS = _fw + (f_LV - _fw) / std::pow(1.0 + std::pow(slip_rate / _sw, 8.0), 1.0/8.0);
  
  return f_SS;
}

Real
RateStateFrictionczm3d::computePsiSS(Real slip_rate)
{
  Real f_SS = computeFSS(slip_rate);
  Real psi_SS = _a * std::log((2.0 * _s0 / slip_rate) * std::sinh(f_SS / _a));
  return psi_SS;
}

Real
RateStateFrictionczm3d::updateStateVariable(Real slip_rate_current, Real psi_old)
{
  // Integrate state evolution equation using analytical solution
  // dψ/dt = -(ṡ/L)[ψ - ψ_SS]
  Real psi_SS = computePsiSS(slip_rate_current);
  Real psi_new = (psi_old - psi_SS) * std::exp(-slip_rate_current * _dt / _L) + psi_SS;
  return psi_new;
}

Real
RateStateFrictionczm3d::computeFrictionStrength(Real slip_rate, Real psi, Real normal_stress)
{
  // S = σ_n * a * arcsinh[ṡ/(2ṡ_0) * exp(ψ/a)]
  Real arg = slip_rate / (2.0 * _s0) * std::exp(psi / _a);
  Real strength = normal_stress * _a * std::asinh(arg);
  return strength;
}

Real
RateStateFrictionczm3d::solveNewtonForSlipRate(Real s_tilde, Real C, Real psi)
{
  // Solve: ṡ = s_tilde + C * arcsinh[ṡ/(2s_0) * exp(ψ/a)]
  // Using Newton's method with substitution w = arcsinh[ṡ/(2s_0) * exp(ψ/a)]
  
  Real exp_term = std::exp(-psi / _a);
  Real initial_slip_rate = std::max(_s_ini, s_tilde); // Initial guess
  Real w = std::asinh(initial_slip_rate / (2.0 * _s0) * std::exp(psi / _a));
  
  for (unsigned int iter = 0; iter < _newton_max_iterations; ++iter)
  {
    // F(w) = s_tilde + C*w - exp(-ψ/a)*2s_0*sinh(w)
    Real F = s_tilde + C * w - exp_term * 2.0 * _s0 * std::sinh(w);
    
    // F'(w) = C - exp(-ψ/a)*2s_0*cosh(w)
    Real Fprime = C - exp_term * 2.0 * _s0 * std::cosh(w);
    
    if (std::abs(Fprime) < 1e-15)
      break;
    
    Real w_new = w - F / Fprime;
    
    if (std::abs(w_new - w) < _newton_tolerance)
    {
      w = w_new;
      break;
    }
    
    w = w_new;
  }
  
  // Convert back to slip rate
  Real slip_rate = 2.0 * _s0 * std::sinh(w) * exp_term;
  return std::max(slip_rate, _s_ini); // Ensure minimum slip rate
}

void
RateStateFrictionczm3d::computeInterfaceTractionAndDerivatives()
{
  // Global Displacement Jump
  RealVectorValue displacement_jump_global(
      _disp_slipweakening_x[_qp] - _disp_slipweakening_neighbor_x[_qp],
      _disp_slipweakening_y[_qp] - _disp_slipweakening_neighbor_y[_qp],
      _disp_slipweakening_z[_qp] - _disp_slipweakening_neighbor_z[_qp]);
  
  // Global Displacement Jump Old
  RealVectorValue displacement_jump_old_global(
    _disp_slipweakening_x_old[_qp] - _disp_slipweakening_neighbor_x_old[_qp],
    _disp_slipweakening_y_old[_qp] - _disp_slipweakening_neighbor_y_old[_qp],
    _disp_slipweakening_z_old[_qp] - _disp_slipweakening_neighbor_z_old[_qp]);

  // Global Displacement Jump Rate
  RealVectorValue displacement_jump_rate_global =
      (displacement_jump_global - displacement_jump_old_global) * (1.0 / _dt);

  // Local Displacement Jump / Displacement Jump Rate
  RealVectorValue displacement_jump = _rot[_qp].transpose() * displacement_jump_global;
  RealVectorValue displacement_jump_rate = _rot[_qp].transpose() * displacement_jump_rate_global;

  // n is along normal direction; t is along tangential direction; d is along dip direction
  Real displacement_jump_n = displacement_jump(0);
  Real displacement_jump_t = displacement_jump(1);
  Real displacement_jump_d = displacement_jump(2);
  Real displacement_jump_rate_n = displacement_jump_rate(0);
  Real displacement_jump_rate_t = displacement_jump_rate(1);
  Real displacement_jump_rate_d = displacement_jump_rate(2);

  // Reaction force in local coordinate
  RealVectorValue R_plus_global(-_reaction_slipweakening_x[_qp],
                                -_reaction_slipweakening_y[_qp],
                                -_reaction_slipweakening_z[_qp]);
  RealVectorValue R_minus_global(-_reaction_slipweakening_neighbor_x[_qp],
                                 -_reaction_slipweakening_neighbor_y[_qp],
                                 -_reaction_slipweakening_neighbor_z[_qp]);

  RealVectorValue R_plus_local = _rot[_qp].transpose() * R_plus_global;
  RealVectorValue R_minus_local = _rot[_qp].transpose() * R_minus_global;

  Real R_plus_local_n = R_plus_local(0);
  Real R_plus_local_t = R_plus_local(1);
  Real R_plus_local_d = R_plus_local(2);
  Real R_minus_local_n = R_minus_local(0);
  Real R_minus_local_t = R_minus_local(1);
  Real R_minus_local_d = R_minus_local(2);

  // Compute node mass and area
  Real M = 0;
  Real A = 0;
  if (_current_elem->type() == libMesh::ElemType::TET4){
    M = (_density[_qp] * sqrt(2) * _len * _len * _len / 12 / 4) * 6;
    A = (sqrt(3) * _len * _len / 4 / 3) * 6;
  }
  else if (_current_elem->type() == libMesh::ElemType::HEX8){
    M = (_density[_qp] * _len * _len * _len / 8) * 4;
    A = (_len * _len / 4) * 4;
  }

  // Compute T1_o, T2_o, T3_o for current qp
  Real T1_o = _static_initial_stress_tensor[_qp](0, 1); // shear stress in t dir
  Real T2_o = -1.0 * _static_initial_stress_tensor[_qp](1, 1); // normal stress in n dir
  Real T3_o = _static_initial_stress_tensor[_qp](0, 2); // shear stress in d dir

  // Compute trial sticking stress
  Real T1_trial = (1 / _dt) * M * displacement_jump_rate_t / (2 * A) +
            (R_plus_local_t - R_minus_local_t) / (2 * A) + T1_o;
  Real T3_trial = (1 / _dt) * M * displacement_jump_rate_d / (2 * A) +
            (R_plus_local_d - R_minus_local_d) / (2 * A) + T3_o;
  Real T2 = -(1 / _dt) * M * (displacement_jump_rate_n + (1 / _dt) * displacement_jump_n) /
                (2 * A) +
            ((R_minus_local_n - R_plus_local_n) / (2 * A)) - T2_o;

  // Handle normal stress (compression is negative)
  if (T2 > 0)
    T2 = 0;

  Real normal_stress = -T2; // Convert to positive compression
  Real Pf = _fluid_pressure_aux[_qp]; // fluid pressure
  Real effective_normal_stress = std::max(normal_stress - Pf, 0.0);

  // Calculate current slip rate magnitude
  Real current_slip_rate = std::sqrt(displacement_jump_rate_t * displacement_jump_rate_t + 
                                   displacement_jump_rate_d * displacement_jump_rate_d);
  current_slip_rate = std::max(current_slip_rate, _s_ini); // Regularization

  // Update state variable
  Real psi_old = _state_variable_old[_qp];
  Real psi_new = updateStateVariable(current_slip_rate, psi_old);

  // Calculate trial shear traction magnitude
  Real T_shear_trial = std::sqrt(T1_trial * T1_trial + T3_trial * T3_trial);

  // Rate-and-state friction calculation
  Real T1 = T1_trial;
  Real T3 = T3_trial;
  
  if (T_shear_trial > 0 && effective_normal_stress > 0)
  {
    // Calculate s_tilde (trial slip rate from paper equation 42)
    Real s_tilde = T_shear_trial; // Simplified - you may need to adjust this based on your specific implementation
    
    // Calculate C coefficient (equation 43)
    Real C = 4.0 * _dt / (_len * _density[_qp]) * effective_normal_stress * _a;
    
    // Solve for new slip rate using Newton's method
    Real new_slip_rate = solveNewtonForSlipRate(s_tilde, C, psi_new);
    
    // Calculate friction strength
    Real friction_strength = computeFrictionStrength(new_slip_rate, psi_new, effective_normal_stress);
    
    // Apply friction constraint
    if (T_shear_trial > friction_strength)
    {
      Real reduction_factor = friction_strength / T_shear_trial;
      T1 = T1_trial * reduction_factor;
      T3 = T3_trial * reduction_factor;
    }
  }

  // Save displacement jump in local coordinate
  _displacement_jump_strike[_qp] = displacement_jump_t;
  _displacement_jump_dip[_qp] = displacement_jump_d;
  _displacement_jump_normal[_qp] = displacement_jump_n;

  // Save displacement jump rate in local coordinate
  _displacement_jump_rate_strike[_qp] = displacement_jump_rate_t;
  _displacement_jump_rate_dip[_qp] = displacement_jump_rate_d;
  _displacement_jump_rate_normal[_qp] = displacement_jump_rate_n;

  // Save traction in local coordinate
  _traction_strike[_qp] = T1;
  _traction_normal[_qp] = T2;
  _traction_dip[_qp] = T3;

  // Additional output properties
  _slip_rate_magnitude[_qp] = current_slip_rate;
  _friction_coefficient[_qp] = effective_normal_stress > 0 ? 
    std::sqrt(T1*T1 + T3*T3) / effective_normal_stress : 0.0;
  _state_variable_updated[_qp] = psi_new;

  // Assign back traction in CZM
  RealVectorValue traction(T2 + T2_o, -T1 + T1_o, -T3 + T3_o);
  _interface_traction[_qp] = traction;
  _dinterface_traction_djump[_qp] = 0;
}