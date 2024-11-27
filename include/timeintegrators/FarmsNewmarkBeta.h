//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "TimeIntegrator.h"
#include "MathUtils.h"

/**
 * Newmark-Beta time integration method
 */
class FarmsNewmarkBeta : public TimeIntegrator
{
public:
  static InputParameters validParams();

  FarmsNewmarkBeta(const InputParameters & parameters);
  virtual ~FarmsNewmarkBeta();

  virtual int order() override { return 1; }
  virtual void computeTimeDerivatives() override;
  virtual void computeADTimeDerivatives(ADReal & ad_u_dot,
                                        const dof_id_type & dof,
                                        ADReal & ad_u_dotdot) const override;
  virtual void postResidual(NumericVector<Number> & residual) override;
  virtual bool overridesSolve() const override { return false; }

protected:
  /**
   * Helper function that actually does the math for computing the time derivative
   */
  template <typename T, typename T2, typename T3, typename T4, typename T5>
  void computeTimeDerivativeHelper(T & u_dot,
                                   const T2 & u_old,
                                   const T3 & u_dot_old,
                                   T4 & u_dotdot,
                                   const T5 & u_dotdot_old,
                                   Real flag,
                                   Real factor,
                                   Real threshold) const;

  virtual Real duDotDuCoeff() const override;

  /// Newmark time integration parameter-beta
  Real _beta;

  /// Newmark time integration parameter-gamma
  Real _gamma;

  /// Inactive time steps
  int _inactive_tsteps;

  /// solution vector for \f$ {du^dotdot}\over{du} \f$
  Real & _du_dotdot_du;

  /// Factor by which the velocity is reduced
  Real _factor;

  /// Threshold value for velocity reduction
  Real _threshold;

};

template <typename T, typename T2, typename T3, typename T4, typename T5>
void
FarmsNewmarkBeta::computeTimeDerivativeHelper(
    T & u_dot, 
    const T2 & u_old, 
    const T3 & u_dot_old, 
    T4 & u_dotdot, 
    const T5 & u_dotdot_old, 
    Real flag, 
    Real factor, 
    Real threshold) const
{
  // compute second derivative
  // according to Newmark-Beta method
  // u_dotdot = first_term - second_term - third_term
  //       first_term = (u - u_old) / beta / dt ^ 2
  //      second_term = u_dot_old / beta / dt
  //       third_term = u_dotdot_old * (1 / 2 / beta - 1)
  u_dotdot -= u_old;
  u_dotdot *= 1.0 / _beta / _dt / _dt;
  MathUtils::addScaled(-1.0 / _beta / _dt, u_dot_old, u_dotdot);
  MathUtils::addScaled(-0.5 / _beta + 1.0, u_dotdot_old, u_dotdot);

  // compute first derivative
  // according to Newmark-Beta method
  // u_dot = first_term + second_term + third_term
  //       first_term = u_dot_old
  //      second_term = u_dotdot_old * (1 - gamma) * dt
  //       third_term = u_dotdot * gamma * dt
  u_dot = u_dot_old;
  MathUtils::addScaled((1.0 - _gamma) * _dt, u_dotdot_old, u_dot);
  MathUtils::addScaled(_gamma * _dt, u_dotdot, u_dot);

  // Force reduce the velocity by a factor of its previous value
  // Velocity reduction logic
  if (flag > 0.0)
  {

    // //outputs
    // if (_verbose)
    // {
    //   _console << "Reducing velocity by factor " << factor << " for DOFs with magnitude less than "
    //            << threshold << std::endl;
    // }

    std::vector<dof_id_type> dof_indices;
    std::vector<Number> values;

    // Calculate the first and last valid indices
    dof_id_type first = u_dot.first_local_index();
    dof_id_type local_size = u_dot.local_size();
    dof_id_type last = first + local_size - 1;

    // Loop over local DOFs safely
    for (dof_id_type i = first; i + 1 <= last; i += 2)
    {
      // Get x and y components
      Number vel_x = u_dot(i);    // x velocity
      Number vel_y = u_dot(i+1);  // y velocity
      
      // Check magnitude of velocity vector
      Real vel_mag = std::sqrt(vel_x*vel_x + vel_y*vel_y);
      
      if (MooseUtils::absoluteFuzzyLessEqual(vel_mag, threshold))
      {
        //std::cout<<"Reducing velocity at DOF "<<i<<" by factor "<<factor<<std::endl;
        // Store indices and reduced values
        dof_indices.push_back(i);
        dof_indices.push_back(i+1);
        values.push_back(factor*vel_x);
        values.push_back(factor*vel_y);
      }
    }

    // Update values using set() for each index-value pair
    for (unsigned int i = 0; i < dof_indices.size(); ++i)
    {
      u_dot.set(dof_indices[i], values[i]);
    }
    
    // Close vector after modifications
    // u_dot.close();
    // u_dotdot.close();
  }
}