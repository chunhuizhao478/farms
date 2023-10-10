/*
Implementation of Damage Evolution Forcing Function (F) :

Strong Form: 

d alpha / dt = F
F = (1-B)[Cd I_2 (xi - xi_o) + D grad^2 alpha] if xi >= xi_o
F = (1-B)[C1 exp(alpha/C2) I2 (xi - xi_o) + D grad^2 alpha] if xi <= xi_o

Weak Form:

int( d(alpha)/dt * v ) - int( (1-B) (Cd I2 (xi - xi_o) * v + D d(alpha)/dx dv/dx ) = 0 if xi >= xi_o

int( d(alpha)/dt * v ) - int( (1-B) (C1 exp(alpha/C2) I2 (xi - xi_o) * v + D d(alpha)/dx dv/dx ) = 0 if xi <= xi_o

- 10/5/2023 Chunhui Zhao

Include power-law correction on Cd (function of strain rate e)
if (e < 1e-4){ Cd = 10 };
else{ Cd = pow(10, log10(1+m*log10(e/1e-4)) ) }

*/

#pragma once

#include "Kernel.h"
#include "Material.h"

class DamageVarForcingFuncDev : public Kernel
{
public:
  static InputParameters validParams();

  DamageVarForcingFuncDev(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();

private:
  /// constant parameters
  Real _Cd_min; //minimum Cd value for small strain (e < 1e-4)
  Real _D;
  Real _C1;
  Real _C2;
  Real _xi_0;
  Real _xi_min;
  Real _xi_max;
  /// scaling the power-law
  Real _scale;
  /// power index
  Real _m;
  /// threshold mechanical strain (Cd remains constant below this value)
  Real _mechanical_strain_rate_threshold;
  /// variable parameters
  const VariableValue & _alpha_old;
  const VariableValue & _B_old;
  const VariableValue & _xi_old;
  const VariableValue & _I2_old;
  /// mechanical strain rate
  const VariableValue & _mechanical_strain_rate;
  /// add options
  int _option;
  Real _Cd_constant;
  /// fix diffusion
  Real _shear_modulus_o;
  Real _lambda_o;
  /// add healing option
  bool _healing;
};