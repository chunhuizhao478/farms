//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

/*
Material Description of Slip Weakening Friction 3d
*/

#include "RateStateFrictionczm3d.h"
#include "InterfaceKernel.h"

registerMooseObject("farmsApp", RateStateFrictionczm3d);

InputParameters
RateStateFrictionczm3d::validParams()
{
  InputParameters params = CZMComputeLocalTractionTotalBase::validParams();
  params.addClassDescription("linear slip weakening traction separation law.");
  params.addRequiredParam<Real>("T1_o", "background shear traction in strike dir");
  params.addRequiredParam<Real>("T2_o", "background normal traction");
  params.addRequiredParam<Real>("T3_o", "background shear traction in dip dir");
  params.addRequiredParam<Real>("len","element length");
  params.addRequiredParam<Real>("f_o","rate-and-state friction coefficients");
  params.addParam<Real>("rsf_a", 0.008, "rate-and-state friction coefficient a (constant value, overridden if rsf_a_var provided)");
  params.addRequiredParam<Real>("rsf_b","rate-and-state friction coefficients");
  params.addRequiredParam<Real>("rsf_L","rate-and-state friction coefficients");
  params.addRequiredParam<Real>("delta_o","slip rate parameter");
  params.addParam<Real>("statevar_init", 1.606238999213454e9, "initial value of state variable (constant, overridden if statevar_init_var provided)");
  params.addRequiredParam<Real>("sliprate_strike_init","initial value of strike slip rate");

  // Optional coupled variables for spatially variable RSF parameters
  params.addCoupledVar("rsf_a_var", "Spatially variable RSF 'a' parameter (optional, overrides rsf_a)");
  params.addCoupledVar("statevar_init_var", "Spatially variable initial state variable (optional, overrides statevar_init)");
  params.addRequiredCoupledVar("disp_x", "displacement in x dir");
  params.addRequiredCoupledVar("disp_y", "displacement in y dir");
  params.addRequiredCoupledVar("disp_z", "displacement in z dir");
  params.addRequiredCoupledVar("vel_x","velocity in x dir");
  params.addRequiredCoupledVar("vel_y","velocity in y dir");
  params.addRequiredCoupledVar("vel_z","velocity in z dir");
  params.addRequiredCoupledVar("reaction_x", "reaction in x dir");
  params.addRequiredCoupledVar("reaction_y", "reaction in y dir");
  params.addRequiredCoupledVar("reaction_z", "reaction in z dir");
  params.addRequiredCoupledVar("Ts_perturb","shear stress perturbation in strike dir at time t");
  return params;
}

RateStateFrictionczm3d::RateStateFrictionczm3d(const InputParameters & parameters)
  : CZMComputeLocalTractionTotalBase(parameters),
    _statevar(declareProperty<Real>("statevar")),
    _statevar_old(getMaterialPropertyOldByName<Real>("statevar")),
    _slipratevar(declareProperty<RealVectorValue>("slipratevar")),
    _slipratevar_old(getMaterialPropertyOldByName<RealVectorValue>("slipratevar")),
    _T1_o(getParam<Real>("T1_o")),
    _T2_o(getParam<Real>("T2_o")),
    _T3_o(getParam<Real>("T3_o")),
    _len(getParam<Real>("len")),
    _f_o(getParam<Real>("f_o")),
    _rsf_a(getParam<Real>("rsf_a")),
    _rsf_b(getParam<Real>("rsf_b")),
    _rsf_L(getParam<Real>("rsf_L")),
    _delta_o(getParam<Real>("delta_o")),
    _statevar_init(getParam<Real>("statevar_init")),
    _sliprate_strike_init(getParam<Real>("sliprate_strike_init")),
    _use_coupled_rsf_a(isParamValid("rsf_a_var") && isCoupled("rsf_a_var")),
    _rsf_a_var(_use_coupled_rsf_a ? &coupledValue("rsf_a_var") : nullptr),
    _use_coupled_statevar_init(isParamValid("statevar_init_var") && isCoupled("statevar_init_var")),
    _statevar_init_var(_use_coupled_statevar_init ? &coupledValue("statevar_init_var") : nullptr),
    _density(getMaterialPropertyByName<Real>(_base_name + "density")),
    _rot(getMaterialPropertyByName<RankTwoTensor>(_base_name + "czm_total_rotation")),
    _disp_x(coupledValue("disp_x")),
    _disp_neighbor_x(coupledNeighborValue("disp_x")),
    _disp_y(coupledValue("disp_y")),
    _disp_neighbor_y(coupledNeighborValue("disp_y")),
    _disp_z(coupledValue("disp_z")),
    _disp_neighbor_z(coupledNeighborValue("disp_z")),
    _vel_x(coupledValue("vel_x")),
    _vel_neighbor_x(coupledNeighborValue("vel_x")),
    _vel_y(coupledValue("vel_y")),
    _vel_neighbor_y(coupledNeighborValue("vel_y")),
    _vel_z(coupledValue("vel_z")),
    _vel_neighbor_z(coupledNeighborValue("vel_z")),
    _reaction_x(coupledValue("reaction_x")),
    _reaction_neighbor_x(coupledNeighborValue("reaction_x")),
    _reaction_y(coupledValue("reaction_y")),
    _reaction_neighbor_y(coupledNeighborValue("reaction_y")),
    _reaction_z(coupledValue("reaction_z")),
    _reaction_neighbor_z(coupledNeighborValue("reaction_z")),
    _disp_x_old(coupledValueOld("disp_x")),
    _disp_neighbor_x_old(coupledNeighborValueOld("disp_x")),
    _disp_y_old(coupledValueOld("disp_y")),
    _disp_neighbor_y_old(coupledNeighborValueOld("disp_y")),
    _disp_z_old(coupledValueOld("disp_z")),
    _disp_neighbor_z_old(coupledNeighborValueOld("disp_z")),
    _Ts_perturb(coupledValue("Ts_perturb")),
    _Ts_perturb_old(coupledValueOld("Ts_perturb")),
    _Tn_debug(declareProperty<Real>("Tn_debug")),
    _Tmag_trial_debug(declareProperty<Real>("Tmag_trial_debug")),
    _T_mag_debug(declareProperty<Real>("T_mag_debug")),
    _sliprate_mag_debug(declareProperty<Real>("sliprate_mag_debug")),
    _Z_debug(declareProperty<Real>("Z_debug")),
    _newton_iters_debug(declareProperty<Real>("newton_iters_debug"))
{

  // only works for small strain
  if (hasBlockMaterialProperty<RankTwoTensor>(_base_name + "strain_increment"))
  {
    mooseError("RateState only works for small strain!");
  }
}

void
RateStateFrictionczm3d::initQpStatefulProperties()
{
  // Use coupled variable for initial state variable if provided, otherwise use constant
  if (_use_coupled_statevar_init)
    _statevar[_qp] = (*_statevar_init_var)[_qp];
  else
    _statevar[_qp] = _statevar_init;

  _slipratevar[_qp](0)     = 0;                     //normal
  _slipratevar[_qp](1)     = _sliprate_strike_init; //strike
  _slipratevar[_qp](2)     = 0;                     //dip

}

void
RateStateFrictionczm3d::computeInterfaceTractionAndDerivatives()
{
  // Get local RSF 'a' parameter (spatially variable or constant)
  Real rsf_a_local = _use_coupled_rsf_a ? (*_rsf_a_var)[_qp] : _rsf_a;

  // Global Displacement Jump
  RealVectorValue displacement_jump_global(
      _disp_x[_qp] - _disp_neighbor_x[_qp],
      _disp_y[_qp] - _disp_neighbor_y[_qp],
      _disp_z[_qp] - _disp_neighbor_z[_qp]);

  // Global Displacement Jump Old
  RealVectorValue displacement_jump_old_global(
    _disp_x_old[_qp] - _disp_neighbor_x_old[_qp],
    _disp_y_old[_qp] - _disp_neighbor_y_old[_qp],
    _disp_z_old[_qp] - _disp_neighbor_z_old[_qp]);

  // Global Displacement Jump Rate
  RealVectorValue displacement_jump_rate_global =
      (displacement_jump_global - displacement_jump_old_global) * (1 / _dt);

  // Local Displacement Jump / Displacement Jump Rate
  RealVectorValue displacement_jump = _rot[_qp].transpose() * displacement_jump_global;
  RealVectorValue displacement_jump_rate = _rot[_qp].transpose() * displacement_jump_rate_global;

  // n is along normal direction; t is along tangential direction; d is along dip direction
  Real displacement_jump_n = displacement_jump(0);
  //Real displacement_jump_t = displacement_jump(1);
  //Real displacement_jump_d = displacement_jump(2);
  Real displacement_jump_rate_n = displacement_jump_rate(0);
  Real displacement_jump_rate_t = displacement_jump_rate(1);
  Real displacement_jump_rate_d = displacement_jump_rate(2);

  // Reaction force in local coordinate
  RealVectorValue R_plus_global(-_reaction_x[_qp],
                                -_reaction_y[_qp],
                                -_reaction_z[_qp]);
  RealVectorValue R_minus_global(-_reaction_neighbor_x[_qp],
                                 -_reaction_neighbor_y[_qp],
                                 -_reaction_neighbor_z[_qp]);

  RealVectorValue R_plus_local = _rot[_qp].transpose() * R_plus_global;
  RealVectorValue R_minus_local = _rot[_qp].transpose() * R_minus_global;

  // n is along normal direction; t is along tangential direction; d is along dip direction
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

  // Compute normal sticking stress
  Real T2 = -(1 / _dt) * M * (displacement_jump_rate_n + (1 / _dt) * displacement_jump_n) /
                (2 * A) +
            ((R_minus_local_n - R_plus_local_n) / (2 * A)) - _T2_o;

  // Compute fault traction
  if (T2 < 0)
  {
  }
  else
  {
    T2 = 0;
  }

  ///Make Tn positive
  Real Tn = abs(T2);

  //*Compute Trial Shear Traction Along Strike Direction at Current Time Step*
  Real Ts_trial = (1 / _dt) * M * displacement_jump_rate_t / (2 * A) +
                  (R_plus_local_t - R_minus_local_t) / (2 * A) + _T1_o + _Ts_perturb[_qp];
  Real Td_trial = (1 / _dt) * M * displacement_jump_rate_d / (2 * A) +
                  (R_plus_local_d - R_minus_local_d) / (2 * A) + _T3_o;
  Real Tmag_trial = sqrt(Ts_trial*Ts_trial+Td_trial*Td_trial);

  //*Compute slip rate magnitude old
  // _slipratevar (n,t,d)
  Real sliprate_mag_old = sqrt(_slipratevar_old[_qp](1)*_slipratevar_old[_qp](1)+_slipratevar_old[_qp](2)*_slipratevar_old[_qp](2));

  //const
  Real c = A * _dt * ( M + M ) / (M * M);

  // FIX: Guard against non-positive state variable in log()
  Real statevar_old_safe = std::max(_statevar_old[_qp], 1e-30);
  Real Z = 0.5 / _delta_o * exp((_f_o + _rsf_b * log(_delta_o * statevar_old_safe/_rsf_L))/rsf_a_local);

  // FIX 1: Guard against Tmag_trial ~ 0
  // When trial traction magnitude is nearly zero, the fault is essentially
  // unloaded — use sticking behavior (keep background traction, no slip change).
  const Real Tmag_trial_tol = 1.0; // 1 Pa tolerance
  if (Tmag_trial < Tmag_trial_tol)
  {
    // No meaningful shear load: stick with previous state
    _interface_traction[_qp] = RealVectorValue(T2 + _T2_o, _T1_o, _T3_o);
    _dinterface_traction_djump[_qp] = 0;
    _statevar[_qp] = _statevar_old[_qp];
    _slipratevar[_qp](1) = 0.0;
    _slipratevar[_qp](2) = 0.0;

    // Debug output
    _Tn_debug[_qp] = Tn;
    _Tmag_trial_debug[_qp] = Tmag_trial;
    _T_mag_debug[_qp] = 0.0;
    _sliprate_mag_debug[_qp] = 0.0;
    _Z_debug[_qp] = Z;
    _newton_iters_debug[_qp] = 0.0;
    return;
  }

  // Compute trial shear direction (safe to divide now)
  Real dir_s = Ts_trial / Tmag_trial;
  Real dir_d = Td_trial / Tmag_trial;

  // =========================================================================
  // Solve for slip rate magnitude using bisection.
  //
  // Residual: f(V) = V + c*Tn*a*asinh(0.5*(V+V_old)*Z) - c*Tmag_trial
  // f is monotonically increasing (f' > 0 always), so bisection is guaranteed
  // to converge. Newton fails when Z is very large (>~1e30) because the
  // asinh flattens to log, making the Jacobian ~1 and Newton oscillates.
  // =========================================================================

  // Lambda for the residual function
  auto resid_func = [&](Real V) -> Real {
    return V + c * Tn * rsf_a_local * asinh(0.5 * (V + sliprate_mag_old) * Z) - c * Tmag_trial;
  };

  // Step 1: Find brackets [V_lo, V_hi] such that f(V_lo) < 0 and f(V_hi) > 0
  Real V_lo, V_hi;
  Real f_at_guess = resid_func(sliprate_mag_old);

  if (abs(f_at_guess) < 1e-14)
  {
    // Already converged — slip rate unchanged
    V_lo = sliprate_mag_old;
    V_hi = sliprate_mag_old;
  }
  else if (f_at_guess > 0)
  {
    // Root is to the left of sliprate_mag_old
    V_hi = sliprate_mag_old;
    Real step = std::max(abs(sliprate_mag_old), 1e-15);
    V_lo = V_hi - step;
    while (resid_func(V_lo) > 0)
    {
      V_hi = V_lo;
      step *= 2.0;
      V_lo = V_hi - step;
    }
  }
  else
  {
    // Root is to the right of sliprate_mag_old
    V_lo = sliprate_mag_old;
    Real step = std::max(abs(sliprate_mag_old), 1e-15);
    V_hi = V_lo + step;
    while (resid_func(V_hi) < 0)
    {
      V_lo = V_hi;
      step *= 2.0;
      V_hi = V_lo + step;
    }
  }

  // Step 2: Bisection to find the root
  int iterr = 0;
  const int max_iter = 200; // bisection converges in ~60 iterations for double precision
  const Real bisect_tol = 1e-12;
  Real solution = 0.5 * (V_lo + V_hi);

  while (iterr < max_iter)
  {
    solution = 0.5 * (V_lo + V_hi);
    Real f_mid = resid_func(solution);

    // Check convergence: relative interval width
    Real interval = V_hi - V_lo;
    if (interval / std::max(abs(solution), 1e-30) < bisect_tol || abs(f_mid) < 1e-14)
      break;

    if (f_mid > 0)
      V_hi = solution;
    else
      V_lo = solution;

    iterr++;
  }

  // Clamp: slip rate magnitude >= 0
  Real sliprate_mag = std::max(solution, 0.0);

  // FIX 3: Floor slip rate for state variable update to prevent L/V -> Inf
  const Real V_floor = 1e-20;
  Real sliprate_for_statevar = std::max(sliprate_mag, V_floor);

  //update state variable
  // we apply separation of variables on first-order linear ODE of statevar
  Real coeff_LoverV = _rsf_L / sliprate_for_statevar;
  Real coeff_exponent = exp(-sliprate_for_statevar*_dt/_rsf_L);
  Real statevar_tplusdt = coeff_LoverV + (_statevar_old[_qp] - coeff_LoverV) * coeff_exponent;

  //*Compute shear traction at time t*
  Real T_mag = Tn * rsf_a_local * asinh( 0.5*(sliprate_mag_old+sliprate_mag) * Z );

  ///Get Components
  Real T1 = T_mag * dir_s;
  Real T3 = T_mag * dir_d;

  // Assign back traction in CZM
  RealVectorValue traction(T2 + _T2_o, -T1 + _T1_o, -T3 + _T3_o);
  _interface_traction[_qp] = traction;
  _dinterface_traction_djump[_qp] = 0;

  // Update statevar, slipratevar (strike, dip)
  _statevar[_qp] = statevar_tplusdt;
  _slipratevar[_qp](1) = sliprate_mag * dir_s;
  _slipratevar[_qp](2) = sliprate_mag * dir_d;

  // Store debug output
  _Tn_debug[_qp] = Tn;
  _Tmag_trial_debug[_qp] = Tmag_trial;
  _T_mag_debug[_qp] = T_mag;
  _sliprate_mag_debug[_qp] = sliprate_mag;
  _Z_debug[_qp] = Z;
  _newton_iters_debug[_qp] = static_cast<Real>(iterr);
}
