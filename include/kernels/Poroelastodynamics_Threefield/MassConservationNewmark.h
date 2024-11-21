#pragma once

#include "Kernel.h"
//Forward Declarations
class MassConservationNewmark : public Kernel
{
public:
  static InputParameters validParams();

  MassConservationNewmark(const InputParameters & parameters);
protected:
virtual Real computeQpResidual();
virtual Real computeQpJacobian();
virtual Real computeQpOffDiagJacobian(unsigned int jvar);
private:
unsigned int _ndisp; // number of displacement components (1D, 2D or 3D)
unsigned int _ux_var; // id of displacement components
unsigned int _uy_var;
unsigned int _uz_var;
const VariableGradient &_grad_ux; // gradient of displacement
const VariableGradient &_grad_uy;
const VariableGradient &_grad_uz;
const VariableGradient &_grad_ux_old; // gradient of previous displacement
const VariableGradient &_grad_uy_old;
const VariableGradient &_grad_uz_old;
const VariableGradient &_grad_vx_old; // gradient of previous velocity
const VariableGradient &_grad_vy_old;
const VariableGradient &_grad_vz_old;
const VariableGradient & _grad_ax_old; // gradient of previous acceleration
const VariableGradient & _grad_ay_old;
const VariableGradient & _grad_az_old;
const Real _beta;
const Real _gamma;
const MaterialProperty<Real> & _coefficient;  /// Biot coefficient
};