//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "CZMComputeLocalTractionTotalBase.h"

/**
 * Implementation of the rate-and-state friction model for 2D problems.
 **/
class RateStateFriction2d : public CZMComputeLocalTractionTotalBase
{
public:
  static InputParameters validParams();
  RateStateFriction2d(const InputParameters & parameters);

protected:
  /// method computing the total traction and its derivatives
  void computeInterfaceTractionAndDerivatives() override;

  /// Material properties
  ///@brief state variable
  MaterialProperty<Real> & _state_variable;
  const MaterialProperty<Real> & _state_variable_old;

  ///@brief slip rate
  MaterialProperty<RealVectorValue> & _slip_rate;
  const MaterialProperty<RealVectorValue> & _slip_rate_old;

  ///@brief slip rate magnitude
  MaterialProperty<Real> & _slip_rate_magnitude;
  const MaterialProperty<Real> & _slip_rate_magnitude_old;

  ///@brief old interface displacement jump
  const MaterialProperty<RealVectorValue> & _interface_displacement_jump_old;

  ///Parameters
  ///@brief friction parameters
  const Real _fo;
  const Real _a;
  const Real _b;
  
  ///@brief reference slip rate V_o
  const Real _slip_rate_ref;
  ///@brief reference state variable theta_o
  const Real _state_variable_ref;
  ///@brief reference length scale L
  const Real _length_scale_ref;
  ///@brief background shear traction T1_o
  const Real _T1_o;
  ///@brief background normal traction T2_o
  const Real _T2_o;

};