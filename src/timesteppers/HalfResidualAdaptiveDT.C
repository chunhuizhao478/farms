#include "HalfResidualAdaptiveDT.h"
#include "FEProblemBase.h"
#include "NonlinearSystemBase.h"
#include "Transient.h"
#include "libmesh/numeric_vector.h"

registerMooseObject("farmsApp", HalfResidualAdaptiveDT);

InputParameters
HalfResidualAdaptiveDT::validParams()
{
  InputParameters params = TimeStepper::validParams();
  params.addClassDescription(
      "Adaptive dt based on ABAQUS-style half-increment residual using a Newmark mid-state test.");

  params.addParam<Real>("initial_dt", 1e-3, "Initial time step size.");
  params.addParam<Real>("half_abs_tol", 1e-8, "Absolute tolerance for ||R_half||.");
  params.addParam<Real>("half_rel_tol", 0.2,  "Relative tolerance factor times ||R_end||.");
  params.addParam<Real>("growth_factor", 1.25, "Multiplier for dt when R_half is acceptable.");
  params.addParam<Real>("shrink_factor", 0.5,  "Multiplier for dt when R_half is large.");
  params.addParam<Real>("max_increase", 2.0,   "Cap on per-step dt growth (ratio).");
  params.addParam<Real>("min_decrease", 0.25,  "Floor on per-step dt shrink (ratio).");
  params.addParam<Real>("dt_min", 0.0,         "Minimum allowed dt.");
  params.addParam<Real>("dt_max", std::numeric_limits<Real>::max(), "Maximum allowed dt.");
  params.addParam<Real>("cutback_factor_at_failure", 0.5, "Factor when step fails to converge.");
  params.addParam<bool>("use_half_residual", true, "If false, skip half residual test (growth/shrink still applied heuristically).");
  params.addParam<unsigned int>("start_adaptive_after", 0, "Enable adaptive (half residual + dt growth/shrink) only after this many completed time steps.");
  return params;
}

HalfResidualAdaptiveDT::HalfResidualAdaptiveDT(const InputParameters & p)
  : TimeStepper(p),
    _half_abs_tol(getParam<Real>("half_abs_tol")),
    _half_rel_tol(getParam<Real>("half_rel_tol")),
    _growth_factor(getParam<Real>("growth_factor")),
    _shrink_factor(getParam<Real>("shrink_factor")),
    _max_increase(getParam<Real>("max_increase")),
    _min_decrease(getParam<Real>("min_decrease")),
    _dt_min(getParam<Real>("dt_min")),
    _dt_max(getParam<Real>("dt_max")),
    _cutback_on_fail(getParam<Real>("cutback_factor_at_failure")),
    _initial_dt(getParam<Real>("initial_dt")),
    _use_half_residual(getParam<bool>("use_half_residual")),
  _start_adaptive_after(getParam<unsigned int>("start_adaptive_after")),
    _next_dt(_initial_dt),
    _first(true),
    _last_end_res_norm(0.0)
{
  if (std::getenv("HALF_RESIDUAL_DT_DEBUG") && _app.processor_id() == 0)
    _console << "[HalfResidualAdaptiveDT] constructed initial_dt=" << _initial_dt << std::endl;
}

Real HalfResidualAdaptiveDT::computeInitialDT()
{
  _first = false;
  _next_dt = std::min(std::max(_initial_dt, _dt_min), _dt_max);
  return _next_dt;
}

Real HalfResidualAdaptiveDT::computeDT()
{
  return std::min(std::max(_next_dt, _dt_min), _dt_max);
}

Real HalfResidualAdaptiveDT::computeFailedDT()
{
  const Real dt_now = getCurrentDT();
  _next_dt = std::min(std::max(dt_now * _cutback_on_fail, _dt_min), _dt_max);
  return _next_dt;
}

void HalfResidualAdaptiveDT::postSolve()
{
  if (!converged())
    return;

  // Early guard: if this is a first-order integrator (no velocity/accel vectors) running in parallel,
  // skip all adaptive logic (diagnostic isolation for sub_app crash).
  {
    auto & nl_guard = _fe_problem.getNonlinearSystemBase(0);
    const bool first_order_guard = (nl_guard.solutionUDot() == nullptr && nl_guard.solutionUDotDot() == nullptr);
    if (first_order_guard && _app.n_processors() > 1)
    {
      if (_app.processor_id() == 0)
        _console << "[HalfResidualAdaptiveDT] postSolve early skip (first-order + parallel)" << std::endl;
      // Conservative mild growth to avoid stagnation
      const Real dt_now_guard = getCurrentDT();
      const Real proposed_guard = std::min(dt_now_guard * _growth_factor, dt_now_guard * _max_increase);
      _next_dt = std::min(proposed_guard, _dt_max);
      return;
    }
  }

  const Real dt_now = getCurrentDT();
  if (dt_now <= 0.0)
    return;

  const bool debug = std::getenv("HALF_RESIDUAL_DT_DEBUG");
  if (debug && _app.processor_id() == 0)
    _console << "[HalfResidualAdaptiveDT] postSolve enter time=" << _time
             << " dt_now=" << dt_now << " step=" << _fe_problem.timeStep() << std::endl;

  // Defer adaptation until after N completed steps (timeStep() counts from 1 after first advance)
  const bool adaptation_enabled = (_fe_problem.timeStep() > _start_adaptive_after);
  if (!adaptation_enabled)
  {
    if (debug && _app.processor_id() == 0)
      _console << "[HalfResidualAdaptiveDT] adaptation deferred (timeStep=" << _fe_problem.timeStep()
               << " <= start_adaptive_after=" << _start_adaptive_after << ")" << std::endl;
    // Just keep current dt (bounded) or apply mild growth if desired (here we hold constant)
    _next_dt = std::min(std::max(dt_now, _dt_min), _dt_max);
    return;
  }

  if (debug && _app.processor_id() == 0)
    _console << "[HalfResidualAdaptiveDT] postSolve computing end residual" << std::endl;
  _last_end_res_norm = computeEndResidualNorm();
  if (debug && _app.processor_id() == 0)
    _console << "[HalfResidualAdaptiveDT] postSolve end residual=" << _last_end_res_norm << std::endl;
  Real Rhalf = 0.0;
  bool half_used = false;
  const bool env_skip = std::getenv("HALF_RESIDUAL_DT_SKIP");
  if (_use_half_residual && !env_skip)
  {
    half_used = true;
    Rhalf = computeHalfResidualNorm(dt_now);
  }

  const Real tol = std::max(_half_abs_tol, _half_rel_tol * _last_end_res_norm);
  const bool ok = (!half_used) ? true : (Rhalf <= tol);

  if (ok)
  {
    const Real proposed = std::min(dt_now * _growth_factor, dt_now * _max_increase);
    _next_dt = std::min(proposed, _dt_max);
  }
  else
  {
    const Real proposed = std::max(dt_now * _shrink_factor, dt_now * _min_decrease);
    _next_dt = std::max(proposed, _dt_min);
  }

  if (debug && _app.processor_id() == 0)
    _console << "[HalfResidualAdaptiveDT] postSolve exit step=" << _fe_problem.timeStep()
             << " end_res=" << _last_end_res_norm
             << (half_used ? (" R_half=" + Moose::stringify(Rhalf)) : " R_half=SKIPPED")
             << " tol=" << tol << " next_dt=" << _next_dt << std::endl;
}

Real HalfResidualAdaptiveDT::computeEndResidualNorm() const
{
  auto & nl = _fe_problem.getNonlinearSystemBase(0);
  // Use existing residual norm from converged solve if available; otherwise recompute
  auto & R = nl.residualGhosted();
  const bool debug = std::getenv("HALF_RESIDUAL_DT_DEBUG");
  const bool skip_end = std::getenv("HALF_RESIDUAL_DT_SKIP_END_RESIDUAL");
  const bool first_order = (nl.solutionUDot() == nullptr && nl.solutionUDotDot() == nullptr);

  // Grab some contextual info early (avoid accessing after potential fault)
  R.close();
  Real current_l2 = R.l2_norm();
  if (debug && _app.processor_id() == 0)
  {
    _console << "[HalfResidualAdaptiveDT] computeEndResidualNorm enter: size=" << R.size()
             << " current_l2=" << current_l2
             << " first_order=" << first_order
             << " procs=" << _app.n_processors()
             << " time(fe_problem)=" << _fe_problem.time()
             << " time(stepper._time)=" << _time
             << std::endl;
  }

  // If vector already non-zero just return it (avoid extra assembly cost)
  if (current_l2 != 0.0)
    return current_l2;

  // Defensive: optionally skip recompute entirely if requested or risky configuration
  if (skip_end)
  {
    if (debug && _app.processor_id() == 0)
      _console << "[HalfResidualAdaptiveDT] computeEndResidualNorm skipping recompute due to env HALF_RESIDUAL_DT_SKIP_END_RESIDUAL" << std::endl;
    return current_l2; // zero
  }

  // In parallel we avoid recomputing end residual to be safe
  if (_app.n_processors() > 1)
  {
    if (debug && _app.processor_id() == 0)
      _console << "[HalfResidualAdaptiveDT] computeEndResidualNorm early return (parallel, l2==0)" << std::endl;
    return current_l2; // zero
  }

  if (debug)
  {
    // Print on every rank (use libMesh::out) to detect which rank dies
    libMesh::out << "[HalfResidualAdaptiveDT][rank=" << _app.processor_id()
                 << "] recompute start" << std::endl;
  }

  // Zero (should already be zero) then recompute residual assembly
  R.zero();
  // Flush changes for safety (some vector impls need close() before assembly add)
  R.close();

  if (debug)
    libMesh::out << "[HalfResidualAdaptiveDT][rank=" << _app.processor_id() << "] about to computeResidual()" << std::endl;

  // Defensive: ensure solution vector size matches residual size before assembly
  if (nl.solution().size() != R.size())
  {
    if (_app.processor_id() == 0)
      _console << "[HalfResidualAdaptiveDT] size mismatch sol=" << nl.solution().size()
               << " R=" << R.size() << " aborting recompute" << std::endl;
    return 0.0;
  }

  // Attempt recompute inside try/catch to emit partial debug if exception (won't catch segfault)
  _fe_problem.computeResidual(nl.solution(), R, nl.number());
  R.close();
  current_l2 = R.l2_norm();

  if (debug)
  {
    libMesh::out << "[HalfResidualAdaptiveDT][rank=" << _app.processor_id()
                 << "] recompute done l2=" << current_l2 << std::endl;
  }

  return current_l2;
}

void HalfResidualAdaptiveDT::buildNewmarkMidState(const NonlinearSystemBase & nl,
                                                  Real dt,
                                                  std::unique_ptr<NumericVector<Number>> & u_mid,
                                                  std::unique_ptr<NumericVector<Number>> & v_mid,
                                                  std::unique_ptr<NumericVector<Number>> & a_mid) const
{
  // Grab endpoints
  const auto & u_n1 = nl.solution();       // at t_{n+1}
  const auto & u_n  = nl.solutionOld();    // at t_n

  const auto * v_n1_p = nl.solutionUDot();
  const auto * v_n_p  = nl.solutionUDotOld();

  const auto * a_n1_p = nl.solutionUDotDot();
  const auto * a_n_p  = nl.solutionUDotDotOld();

  // If velocities/accelerations are unavailable (e.g. first-order integrator), fall back to simple midpoint
  if (!(v_n1_p && v_n_p && a_n1_p && a_n_p))
  {
    u_mid = u_n.clone();
    u_mid->zero();
    // arithmetic midpoint in displacement space
    u_mid->add(0.5, u_n);
    u_mid->add(0.5, u_n1);
    v_mid.reset();
    a_mid.reset();
    return;
  }

  const auto & v_n1 = *v_n1_p;
  const auto & v_n  = *v_n_p;
  const auto & a_n1 = *a_n1_p;
  const auto & a_n  = *a_n_p;

  // a_{1/2} = 0.5*(a_n + a_{n+1})
  a_mid = a_n.clone();
  a_mid->zero();
  a_mid->add(0.5, a_n);
  a_mid->add(0.5, a_n1);

  // v_{1/2}
  v_mid = v_n.clone();
  v_mid->zero();
  v_mid->add(0.5, v_n);
  v_mid->add(0.5, v_n1);

  // u_{1/2}
  u_mid = u_n.clone();
  u_mid->zero();
  u_mid->add(1.0, u_n);                         // u_n
  if (dt > 0.0)
  {
    std::unique_ptr<NumericVector<Number>> tmp = v_n.clone();
    u_mid->add(0.5 * dt, v_n);                  // 0.5*dt*v_n
    tmp->zero();
    tmp->add(1.0, a_n);
    tmp->add(1.0, a_n1);
    u_mid->add(0.125 * dt * dt, *tmp);          // (dt^2/8)*(a_n+a_{n+1})
  }
}

Real HalfResidualAdaptiveDT::computeHalfResidualNorm(Real dt)
{
  // If no velocity/accel (first-order integrator) AND running in parallel, skip for now (stability)
  auto & nl = _fe_problem.getNonlinearSystemBase(0);
  if (!_use_half_residual)
    return 0.0;
  if (_app.n_processors() > 1 && std::getenv("HALF_RESIDUAL_DT_DISABLE_PARALLEL"))
  {
    if (std::getenv("HALF_RESIDUAL_DT_DEBUG") && _app.processor_id() == 0)
      _console << "[HalfResidualAdaptiveDT] skipping half residual (env disable in parallel)" << std::endl;
    return 0.0;
  }
  if (nl.solutionUDot() == nullptr && nl.solutionUDotOld() == nullptr && _app.n_processors() > 1)
  {
    if (std::getenv("HALF_RESIDUAL_DT_DEBUG") && _app.processor_id() == 0)
      _console << "[HalfResidualAdaptiveDT] skipping half residual (first-order + parallel)" << std::endl;
    return 0.0;
  }
  // Skip if no actual time advance (initial strange step where time==time_old)
  if (std::abs(_time - _time_old) < 1e-20)
  {
    if (std::getenv("HALF_RESIDUAL_DT_DEBUG") && _app.processor_id() == 0)
      _console << "[HalfResidualAdaptiveDT] skipping half residual (time_old==time)" << std::endl;
    return 0.0;
  }
  if (std::getenv("HALF_RESIDUAL_DT_DEBUG"))
    _console << "[HalfResidualAdaptiveDT] computeHalfResidualNorm start dt=" << dt
             << " time_old=" << _time_old << " time_end=" << _time << std::endl;
  
  std::unique_ptr<NumericVector<Number>> u_mid, v_mid, a_mid;
  if (std::getenv("HALF_RESIDUAL_DT_DEBUG"))
    libMesh::out << "[HalfResidualAdaptiveDT][rank=" << _app.processor_id() << "] midstate build begin" << std::endl;
  buildNewmarkMidState(nl, dt, u_mid, v_mid, a_mid);
  if (std::getenv("HALF_RESIDUAL_DT_DEBUG"))
    libMesh::out << "[HalfResidualAdaptiveDT][rank=" << _app.processor_id() << "] midstate build done" << std::endl;

  // Compute t_mid for logging only; avoid mutating global time to keep parallel/MultiApp safe
  const Real t_mid  = 0.5 * (_time_old + _time);
  if (std::getenv("HALF_RESIDUAL_DT_DEBUG"))
    libMesh::out << "[HalfResidualAdaptiveDT][rank=" << _app.processor_id() << "] eval t_mid=" << t_mid << std::endl;

  // Assemble residual at t_{n+1/2} using the mid-state vector directly, avoiding in-place overwrites
  auto & R = nl.residualGhosted();
  if (std::getenv("HALF_RESIDUAL_DT_DEBUG"))
    libMesh::out << "[HalfResidualAdaptiveDT][rank=" << _app.processor_id() << "] residual assembly begin" << std::endl;
  R.zero();
  // Use u_mid as the trial solution; v_mid/a_mid are used only to form u_mid in this approach
  _fe_problem.computeResidual(*u_mid, R, nl.number());
  R.close();
  const Real Rnorm = R.l2_norm();
  if (std::getenv("HALF_RESIDUAL_DT_DEBUG"))
    libMesh::out << "[HalfResidualAdaptiveDT][rank=" << _app.processor_id() << "] residual assembly done R_half=" << Rnorm << std::endl;

  // Optional debug print (enable by setting env HALF_RESIDUAL_DT_DEBUG)
  if (std::getenv("HALF_RESIDUAL_DT_DEBUG"))
  {
    _console << "[HalfResidualAdaptiveDT] computeHalfResidualNorm end dt=" << dt
              << " R_half=" << Rnorm << " t_mid=" << t_mid << std::endl;
  }

  return Rnorm;
}
