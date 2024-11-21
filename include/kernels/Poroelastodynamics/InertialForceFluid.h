
#pragma once

#include "TimeKernel.h"
#include "ADTimeKernel.h"
#include "Material.h"

// Forward Declarations
class TimeIntegrator;

class InertialForceFluid : public TimeKernel
{
public:
  static InputParameters validParams();

  InertialForceFluid(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  const VariableValue * _u_old;
  const VariableValue * _vel_old;
  const VariableValue * _accel_old;
  const bool _has_beta;
  const bool _has_gamma;
  const Real _beta;
  const Real _gamma;
  const bool _has_velocity;
  const bool _has_acceleration;
  const Real _eta;
  const Real _alpha;

  // Velocity and acceleration calculated by time integrator
  const VariableValue * _u_dot_factor_dof;
  const VariableValue * _u_dotdot_factor_dof;
  const VariableValue * _u_dot_factor;
  const VariableValue * _u_dotdot_factor;
  const VariableValue * _u_dot_old;
  const VariableValue * _du_dot_du;
  const VariableValue * _du_dotdot_du;

  /// The fluid component index
  const unsigned int _fluid_component;

  /// PorousFlowDictator UserObject
  const PorousFlowDictator & _dictator;

  /// Number of fluid phases
  const unsigned int _num_phases;

  /// Fluid density for each phase (at the qp)
  const MaterialProperty<std::vector<Real>> & _fluid_density_qp;

  /// Porosity at the  (at the qp)
  const MaterialProperty<Real> & _porosity;
  
  /// The TimeIntegrator
  TimeIntegrator & _time_integrator;
 
};