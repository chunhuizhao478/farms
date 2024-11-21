#pragma once

#include "TimeKernel.h"

// Forward Declarations
class TimeIntegrator;
//Forward Declarations
class INSmassSolid : public TimeKernel
{
public:
  static InputParameters validParams();

  INSmassSolid(const InputParameters & parameters);

protected:
virtual Real computeQpResidual() override;
virtual Real computeQpJacobian() override;
virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;
private:
const VariableGradient * _grad_ux_dot;
const VariableGradient * _grad_uy_dot;
const VariableGradient * _grad_uz_dot;
const MaterialProperty<Real> & _coefficient;  /// Biot coefficient
unsigned int _ux_var; // id of displacement components
unsigned int _uy_var;
unsigned int _uz_var;
const VariableValue * _d_grad_ux_dot_dv;
const VariableValue * _d_grad_uy_dot_dv;
const VariableValue * _d_grad_uz_dot_dv;
MooseVariable & _vx_var;
MooseVariable & _vy_var;
MooseVariable & _vz_var;
const VariablePhiGradient & _vectorx_phi;
const VariablePhiGradient & _vectory_phi;
const VariablePhiGradient & _vectorz_phi;
/// The TimeIntegrator
 TimeIntegrator & _time_integrator;
};