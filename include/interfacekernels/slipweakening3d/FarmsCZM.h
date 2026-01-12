//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "InterfaceKernel.h"

/**
 * Implement Linear Slip Weakening Friction Law Interface Kernel
 */
class FarmsCZM : public InterfaceKernel
{
public:
  static InputParameters validParams();
  FarmsCZM(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual(Moose::DGResidualType type) override;
  virtual Real computeQpJacobian(Moose::DGJacobianType type) override;

  const MaterialProperty<RealVectorValue> & _traction_on_interface;
  const MaterialProperty<RealTensorValue> & _material_tangent_modulus_on_interface;

  /**
   * the displacement component this kernel is operating on (0=x, 1=y, 2 =z)
   */
  unsigned int _component;
  
  /**
   * number of displacement components
   */
  const unsigned int _ndisp;

  std::vector<unsigned int> _disp_var;
};
