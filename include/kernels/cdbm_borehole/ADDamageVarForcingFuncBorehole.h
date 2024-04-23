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

#include "ADKernel.h"
#include "Material.h"

class ADDamageVarForcingFuncBorehole : public ADKernel
{
public:
  static InputParameters validParams();

  ADDamageVarForcingFuncBorehole(const InputParameters & parameters);

protected:

  virtual ADReal computeQpResidual();

private:
  /// constant parameters
  ADReal _Cd;
  ADReal _D;
  ADReal _C1;
  ADReal _C2;
  ADReal _xi_0;
  ADReal _xi_min;
  ADReal _xi_max;
  /// variable parameters
  const ADVariableValue & _alpha_old;
  const ADVariableValue & _B_old;
  const ADVariableValue & _xi_old;
  const ADVariableValue & _I2_old;
};
