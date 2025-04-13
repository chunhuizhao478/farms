//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "AuxKernel.h"

/**
 * postprocess interface material property, this auxkernel may need to define only on interfaces (boundaries) 
 * Created By Chunhui Zhao, Apr9, 2025
 */
class FarmsMaterialRealAux : public AuxKernel
{
public:
  static InputParameters validParams();

  FarmsMaterialRealAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  /// Request material name
  std::string _material_property_name;

  /// The value being set for the current node/element
  const MaterialProperty<Real> & _displacement_jump_global_x;
  const MaterialProperty<Real> & _displacement_jump_global_y;
  const MaterialProperty<Real> & _traction_global_x; 
  const MaterialProperty<Real> & _traction_global_y;
  const MaterialProperty<Real> & _displacement_jump_global_x_old;
  const MaterialProperty<Real> & _displacement_jump_global_y_old;
  const MaterialProperty<Real> * _fractal_shear_stress;
 
  const VariableValue & _ini_shear_sts;
  const VariableValue & _ini_normal_sts;

  /// normals at quadrature points
  const MooseArray<Point> & _normals;

  /// use fractal shear stress
  bool _use_fractal_shear_stress;

};