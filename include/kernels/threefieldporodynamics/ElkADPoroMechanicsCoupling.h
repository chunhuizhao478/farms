//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ADKernel.h"

// Forward Declarations

/**
 * Created by Chunhui Zhao, Jul 5th, 2024
 * ElkADPoroMechanicsCoupling computes -coefficient*porepressure*grad_test[component]
 */
class ElkADPoroMechanicsCoupling : public ADKernel
{
public:
  static InputParameters validParams();

  ElkADPoroMechanicsCoupling(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual();

private:

  /// whether or not multiply biot coefficient by weak form
  const bool _multiply_biot_coefficient;

  /// Biot coefficient
  const ADMaterialProperty<Real> & _coefficient;

  const ADVariableValue & _porepressure;

  unsigned int _porepressure_var_num;

  /// An integer corresponding to the direction this kernel acts in
  unsigned int _component;
};