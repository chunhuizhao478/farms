//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

/*
 * Perfectly Matched Layer (PML) damping kernel for absorbing boundaries
 * Implements damping term for time-domain PML in elastic wave propagation
 * Created for farms_aftershock project
 */

#pragma once

#include "Kernel.h"

/**
 * PMLDamping adds spatially-varying damping in the PML region
 * to absorb outgoing waves without reflection.
 *
 * The damping is applied as: -d(x,y) * rho * du/dt
 * where d(x,y) is the damping coefficient computed by PMLCoefficientMaterial
 */
class PMLDamping : public Kernel
{
public:
  static InputParameters validParams();

  PMLDamping(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

private:
  /// Velocity (time derivative of displacement)
  const VariableValue & _u_dot;

  /// Density material property
  const MaterialProperty<Real> & _density;

  /// PML damping coefficient (spatially varying)
  const MaterialProperty<Real> & _pml_damping_coeff;
};
