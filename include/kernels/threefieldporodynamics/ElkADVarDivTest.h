//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "FEProblemBase.h"
#include "ADKernel.h"

/**
 *  Weak form contribution corresponding to (-p div dw )
 */
class ElkADVarDivTest : public ADKernel
{
public:
  static InputParameters validParams();

  ElkADVarDivTest(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  /// value of the coupled scalar variable
  const ADVariableValue & _p;

  /// Optional coupled concentration variable
  const MooseVariable * _wx;

  /// Optional coupled concentration variable
  const MooseVariable * _wy;

  /// Gradient of the shape function wx
  const VariablePhiGradient & _grad_wx_phi;

  /// Gradient of the shape function wy
  const VariablePhiGradient & _grad_wy_phi;

};