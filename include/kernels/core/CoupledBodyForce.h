//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "GenericKernel.h"

/**
 * CoupledBodyForce kernel similar to BodyForce but scales the body force by a
 * material property (e.g., depth-dependent density). The material property name
 * is provided by input parameter `density_property_name`.
 */
template <bool is_ad>
class CoupledBodyForceTempl : public GenericKernel<is_ad>
{
public:
  static InputParameters validParams();

  CoupledBodyForceTempl(const InputParameters & parameters);

protected:
  virtual GenericReal<is_ad> computeQpResidual() override;

  /// Constant scale factor (optional)
  const Real & _scale;

  /// Optional function value
  const Function & _function;

  /// Optional Postprocessor value
  const PostprocessorValue & _postprocessor;

  /// AD/non-AD version of the quadrature point coordinates
  const MooseArray<Moose::GenericType<Point, is_ad>> * _generic_q_point;

  /// Density-like material property used as multiplicative scale
  const MaterialProperty<Real> & _density;

  usingGenericKernelMembers;
};

typedef CoupledBodyForceTempl<false> CoupledBodyForce;
typedef CoupledBodyForceTempl<true> ADCoupledBodyForce;
