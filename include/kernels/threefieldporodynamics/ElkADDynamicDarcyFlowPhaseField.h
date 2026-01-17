//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

/**
 * ElkADDynamicDarcyFlowPhaseField - Dynamic Darcy flow kernel for phase-field damage coupling
 *
 * Implements: rho^f a^s + rho^f tau_t / phi * a^f + mu_f / kappa * v^f
 *
 * Uses RankTwoTensor for permeability (AD-compatible with phase-field materials)
 *
 * rho^f: fluid density
 * a^s: solid acceleration (Newmark)
 * tau_t: tortosity
 * phi: porosity (damage-dependent)
 * a^f: fluid acceleration (Newmark)
 * mu_f: viscosity
 * kappa: permeability (damage-dependent, RankTwoTensor)
 * v^f: Darcy velocity (Newmark)
 */
#pragma once

#include "ADKernel.h"

class ElkADDynamicDarcyFlowPhaseField : public ADKernel
{
public:
  static InputParameters validParams();
  ElkADDynamicDarcyFlowPhaseField(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

private:
  /// Fluid density
  const ADMaterialProperty<Real> & _rhof;

  /// Porosity (damage-dependent)
  const ADMaterialProperty<Real> & _nf;

  /// Tortosity
  const ADMaterialProperty<Real> & _taut;

  /// Viscosity
  const ADMaterialProperty<Real> & _muf;

  /// Permeability tensor (RankTwoTensor for AD compatibility with phase-field)
  const ADMaterialProperty<RankTwoTensor> & _kappa;

  /// Current skeleton displacement
  const ADVariableValue & _us;

  /// Old skeleton displacement
  const VariableValue & _us_old;

  /// Old skeleton velocity
  const VariableValue & _vs_old;

  /// Old skeleton acceleration
  const VariableValue & _as_old;

  /// Old fluid displacement (primary variable)
  const VariableValue & _fluiddisp_old;

  /// Old fluid velocity
  const VariableValue & _fluidvel_old;

  /// Old fluid acceleration
  const VariableValue & _fluidaccel_old;

  /// Newmark beta parameter
  const Real _beta;

  /// Newmark gamma parameter
  const Real _gamma;

  /// Component direction (0=x, 1=y, 2=z)
  const int _component;
};
