//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "FarmsActuallyExplicitEuler.h"

/**
 * Implements a truly explicit (no nonlinear solve) Central Difference time
 * integration scheme.
 */
class FarmsCentralDifference : public FarmsActuallyExplicitEuler
{
public:
  FarmsCentralDifference(const InputParameters & parameters);

  static InputParameters validParams();

  virtual void initialSetup() override;

  virtual int order() override { return 2; }
  virtual void computeTimeDerivatives() override;
  void computeADTimeDerivatives(ADReal & ad_u_dot,
                                const dof_id_type & dof,
                                ADReal & ad_u_dotdot) const override;

protected:
  virtual Real duDotDuCoeff() const override;

  /// solution vector for \f$ {du^dotdot}\over{du} \f$
  Real & _du_dotdot_du;

  /// The older solution
  const NumericVector<Number> & _solution_older;

  /**
   * Helper function that actually does the math for computing the time derivative
   */
  template <typename T, typename T2, typename T3, typename T4>
  void
  computeTimeDerivativeHelper(T & u_dot, T2 & u_dotdot, const T3 & u_old, const T4 & u_older) const;
};

template <typename T, typename T2, typename T3, typename T4>
void
FarmsCentralDifference::computeTimeDerivativeHelper(T & u_dot,
                                               T2 & u_dotdot,
                                               const T3 & u_old,
                                               const T4 & u_older) const
{
  // computing first derivative
  // using the Central Difference method
  // u_dot_old = (first_term - second_term) / 2 / dt
  //       first_term = u
  //      second_term = u_older
  u_dot -= u_older; // 'older than older' solution
  u_dot *= 1.0 / (2.0 * _dt);

  // computing second derivative
  // using the Central Difference method
  // u_dotdot_old = (first_term - second_term + third_term) / dt / dt
  //       first_term = u
  //      second_term = 2 * u_old
  //       third_term = u_older
  u_dotdot -= u_old;
  u_dotdot -= u_old;
  u_dotdot += u_older;
  u_dotdot *= 1.0 / (_dt * _dt);
}