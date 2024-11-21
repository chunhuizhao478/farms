//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "TimeKernel.h"
#include "Material.h"

// Forward Declarations
class TimeIntegrator;

class FluidStorage : public TimeKernel
{
public:
  static InputParameters validParams();

  FluidStorage(const InputParameters & parameters);

  virtual void computeJacobian() override;

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  bool _lumping;
  const MaterialProperty<Real> & _Biot_modulus;  /// Biot coefficient
  const VariableValue * _u_dot_factor;
  const VariableValue * _u_old;
  const VariableValue * _du_dot_du;
  const VariableValue * _u_dot_factor_dof;
  /// The TimeIntegrator
  TimeIntegrator & _time_integrator;
};