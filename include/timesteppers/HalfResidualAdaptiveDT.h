#pragma once
#include "TimeStepper.h"

class HalfResidualAdaptiveDT : public TimeStepper
{
public:
  static InputParameters validParams();
  HalfResidualAdaptiveDT(const InputParameters & params);

  Real computeInitialDT() override;
  Real computeDT() override;
  Real computeFailedDT() override;
  void postSolve() override;

protected:
  // Core computations
  Real computeHalfResidualNorm(Real dt);
  Real computeEndResidualNorm() const;

  void buildNewmarkMidState(const class NonlinearSystemBase & nl,
                            Real dt,
                            std::unique_ptr<class NumericVector<Number>> & u_mid,
                            std::unique_ptr<class NumericVector<Number>> & v_mid,
                            std::unique_ptr<class NumericVector<Number>> & a_mid) const;

  // Tunables
  const Real _half_abs_tol;
  const Real _half_rel_tol;
  const Real _growth_factor;
  const Real _shrink_factor;
  const Real _max_increase;
  const Real _min_decrease;
  const Real _dt_min;
  const Real _dt_max;
  const Real _cutback_on_fail;
  const Real _initial_dt;
  const bool _use_half_residual;
  const unsigned int _start_adaptive_after; // number of completed time steps before enabling half residual logic

  // State
  Real _next_dt;
  bool _first;
  Real _last_end_res_norm;
};