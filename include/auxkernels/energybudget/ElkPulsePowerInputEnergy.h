//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "AuxKernel.h"

/**
 * Compute Pulse Power Input Energy
 */
class ElkPulsePowerInputEnergy : public AuxKernel
{
public:
  static InputParameters validParams();

  ElkPulsePowerInputEnergy(const InputParameters & parameters);

protected:
  virtual Real computeValue();

  /// option
  int _option;

  /// Function being used to simulate pulse power
  const Function & _func;

  /// Const Value
  Real _confinement_pressure;

  /// number of components in _disp
  unsigned int _ncomp;

  /// displacement variable disp_x, disp_y, disp_z
  std::vector<const VariableValue *> _disp;

  /// old displacement variable disp_x, disp_y, disp_z
  std::vector<const VariableValue *> _disp_old;

  /// normals at quadrature points
  const MooseArray<Point> & _normals;

};