//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "PoroCZMInterfaceKernelbase.h"

/// DG Poro cohesive zone model kernel for the small strain formulation.
/// This kernel assummes the traction sepration law  depends from the
/// displacement jump and fluid velocity jump.
class PoroCZMInterfaceKernelSmallStrain : public PoroCZMInterfaceKernelbase
{
public:
  static InputParameters validParams();
  PoroCZMInterfaceKernelSmallStrain(const InputParameters & parameters);

protected:
  Real computeDResidualDDisplacement(const unsigned int & component_j,
                                     const Moose::DGJacobianType & type) const override;
  // Real computeDResidualDFlux(const unsigned int & component_j,
  //                                    const Moose::DGJacobianType & type) const override;
  // Real computeDResidualDPressure(const Moose::DGJacobianType & type) const override; 
};