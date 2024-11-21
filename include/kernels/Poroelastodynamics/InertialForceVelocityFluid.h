
#pragma once

#include "TimeKernel.h"
#include "ADTimeKernel.h"
#include "Material.h"

// Forward Declarations
class TimeIntegrator;

class InertialForceVelocityFluid : public TimeKernel
{
public:
  static InputParameters validParams();

  InertialForceVelocityFluid(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  // Velocity and acceleration calculated by time integrator
  const VariableValue * _u_dot_factor;
  const VariableValue * _du_dot_du;

  /// The fluid component index
  const unsigned int _fluid_component;

  /// Apparent density
  const Real _apparent_density;

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