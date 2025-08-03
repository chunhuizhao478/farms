#pragma once
//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADKernel.h"

/**
 * This kernel implements a local equivalent strain force using automatic differentiation
 */
class ADCoupledElkLocalEqstrainForce : public ADKernel
{
public:
  static InputParameters validParams();

  ADCoupledElkLocalEqstrainForce(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

private:
  const ADMaterialProperty<Real> & _eqstrain_local;
  const Real _length_scale;
  const Real _kappa_i;
  const Real _c0; // Minimum value of the gradient activity parameter for the equivalent strain
};