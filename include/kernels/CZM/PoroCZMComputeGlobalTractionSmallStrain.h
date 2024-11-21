//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "PoroCZMComputeGlobalTractionBase.h"
/**
 * This class uses the interface traction and its derivatives w.r.t. the interface displacment and fluid velocity jump
 * to compute their respective values in global coordinates under the small strain assumption. The
 * values computed by this object are used by the PoroCZMInterfaceKernelSmallStrain to add the proper
 * residual to the system and to compute the analytic Jacobian.
 */

class PoroCZMComputeGlobalTractionSmallStrain : public PoroCZMComputeGlobalTractionBase
{
public:
  static InputParameters validParams();
  PoroCZMComputeGlobalTractionSmallStrain(const InputParameters & parameters);

protected:
  /// method used to compute the traction and it's derivatives in global coordinates.
  void computeEquilibriumTractionAndPressureAndDerivatives() override;
};