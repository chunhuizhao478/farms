#pragma once

#include "TimeKernel.h"
#include "Material.h"
#include "InertialForce.h"

//Forward Declarations
class TimeIntegrator;
//Forward Declarations
class FluidInertialForce : public TimeKernel
{
public:
  static InputParameters validParams();

  FluidInertialForce(const InputParameters & parameters);

protected:
virtual Real computeQpResidual() override;
virtual Real computeQpJacobian() override;
private:
bool _lumping;
const MaterialProperty<Real> & _rhof; // fluid density
const MaterialProperty<Real> & _nf; // porosity
const VariableValue * _us_dot_dot; // skeleton acceleration for explicit
const VariableValue * _dus_dot_dot_du;
const VariableValue * _u_dot_factor;
const VariableValue * _du_dot_du;
const VariableValue * _u_dot_factor_dof;
/// The TimeIntegrator
const TimeIntegrator & _time_integrator;
};
