//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "PorousFlowPermeabilityBase.h"

/**
 * Material designed to provide a constant permeability tensor
 */
class ElkPorousFlowPermeabilityDamaged : public PorousFlowPermeabilityBase
{
public:
  static InputParameters validParams();

  ElkPorousFlowPermeabilityDamaged(const InputParameters & parameters);

protected:
  void computeQpProperties() override;

  /// Effective permeability tensor
  const MaterialProperty<RealTensorValue> & _effective_perm;

};