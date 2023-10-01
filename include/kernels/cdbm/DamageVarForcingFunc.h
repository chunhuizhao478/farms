/*
Implementation of Damage Evolution Forcing Function (F) :

Strong Form: 

d alpha / dt = F
F = (1-B)[Cd I_2 (xi - xi_o) + D grad^2 alpha] if xi >= xi_o
F = (1-B)[C1 exp(alpha/C2) I2 (xi - xi_o) + D grad^2 alpha] if xi <= xi_o

Weak Form:

int( d(alpha)/dt * v ) - int( (1-B) (Cd I2 (xi - xi_o) * v + D d(alpha)/dx dv/dx ) = 0 if xi >= xi_o

int( d(alpha)/dt * v ) - int( (1-B) (C1 exp(alpha/C2) I2 (xi - xi_o) * v + D d(alpha)/dx dv/dx ) = 0 if xi <= xi_o

*/

#pragma once

#include "Kernel.h"
#include "Material.h"

class DamageVarForcingFunc : public Kernel
{
public:
  static InputParameters validParams();

  DamageVarForcingFunc(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();

private:
  /// constant parameters
  Real _Cd;
  Real _D;
  Real _C1;
  Real _C2;
  Real _xi_0;
  Real _xi_min;
  Real _xi_max;
  /// variable parameters
  const VariableValue & _alpha_old;
  const VariableValue & _B_old;
  const VariableValue & _xi_old;
  const VariableValue & _I2_old;
};
