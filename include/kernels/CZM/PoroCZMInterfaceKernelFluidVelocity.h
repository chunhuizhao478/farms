//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "PoroCZMInterfaceKernelFluxbase.h"

/// DG Poro cohesive zone model kernel for the small strain formulation.
/// This kernel assummes the traction sepration law  depends from the
/// displacement jump and fluid velocity jump.
class PoroCZMInterfaceKernelFluidVelocity : public PoroCZMInterfaceKernelFluxbase
{
public:
  static InputParameters validParams();
  PoroCZMInterfaceKernelFluidVelocity(const InputParameters & parameters);

protected:
  Real computeDResidualDFlux(const unsigned int & component_j,
                                     const Moose::DGJacobianType & type) const override;
};