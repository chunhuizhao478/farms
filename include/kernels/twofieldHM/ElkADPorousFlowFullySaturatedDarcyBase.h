//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

/**
 * Created by Chunhui Zhao, Aug 11th, 2025
 * Implement flux * grad(test)
 * flux: k / mu_f grad(p)
 * k: permeability
 * mu_f: fluid viscosity
 * p: pore pressure
 */
#pragma once

#include "ADKernel.h"
#include "Material.h"

//Forward Declarations

class ElkADPorousFlowFullySaturatedDarcyBase : public ADKernel
{
public:
    static InputParameters validParams();
    ElkADPorousFlowFullySaturatedDarcyBase(const InputParameters & parameters);
    
protected:
    virtual ADReal computeQpResidual() override;

private:
  /// fluid viscosity
  const ADMaterialProperty<Real> & _viscosity;

  /// Material permeability
  const ADMaterialProperty<RealTensorValue> & _permeability;    
};