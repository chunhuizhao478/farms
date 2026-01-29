//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "TimeStepper.h"
#include "PostprocessorInterface.h"

/**
 * SEASAdaptiveDT implements adaptive time stepping for SEAS simulations.
 *
 * The time step is adapted based on the maximum slip rate on the fault:
 * - During aseismic periods (max_V < V_seismic): use larger time steps
 * - During seismic periods (max_V >= V_seismic): use smaller time steps
 *
 * The time step formula is based on the characteristic slip time:
 *   dt = C * Dc / max_V
 *
 * with bounds dt_min <= dt <= dt_max
 *
 * Reference: SCEC SEAS BP2 benchmark specification
 */
class SEASAdaptiveDT : public TimeStepper, public PostprocessorInterface
{
public:
  static InputParameters validParams();

  SEASAdaptiveDT(const InputParameters & parameters);

protected:
  virtual Real computeInitialDT() override;
  virtual Real computeDT() override;

  /// Postprocessor providing maximum slip rate on the fault
  const PostprocessorValue & _max_slip_rate;

  /// Critical slip distance Dc
  const Real _Dc;

  /// Safety factor for time step calculation
  const Real _C;

  /// Minimum allowed time step
  const Real _dt_min;

  /// Maximum allowed time step
  const Real _dt_max;

  /// Threshold slip rate for seismic vs aseismic (default 1e-3 m/s)
  const Real _V_seismic;

  /// Time step during seismic phase
  const Real _dt_seismic;

  /// Initial time step
  const Real _initial_dt;

  /// Growth factor for dt during aseismic phase
  const Real _growth_factor;
};
