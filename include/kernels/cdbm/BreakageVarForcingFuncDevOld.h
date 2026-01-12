/*
Update Breakage Variable Using Runge-Kutta own solver

- 10/5/2023 Chunhui Zhao

Include power-law correction on Cd (function of strain rate e)
if (e < 1e-4){ Cd = 10 };
else{ Cd = pow(10, log10(1+m*log10(e/1e-4)) ) }

Cb = CdCb_multiplier * Cd

*/

#pragma once

#include "Kernel.h"
#include "Material.h"

class BreakageVarForcingFuncDevOld : public Kernel
{
public:
  static InputParameters validParams();

  BreakageVarForcingFuncDevOld(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();

private:
    /// constant parameters
    Real _Cd_min; //minimum Cd value for small strain (e < 1e-4)
    /// multiplier between Cb and Cbh
    Real _CBCBH_multiplier;
    Real _a0;
    Real _a1;
    Real _a2;
    Real _a3;
    Real _xi_min;
    Real _xi_max;
    Real _xi_d;
    Real _xi_0;
    Real _xi_1;
    Real _beta_width;
    /// scaling the power-law
    Real _scale;
    /// power index
    Real _m;
    /// threshold mechanical strain (Cd remains constant below this value)
    Real _mechanical_strain_rate_threshold;
    /// multiplier between Cd and Cb
    Real _CdCb_multiplier;
    /// variable parameters
    const VariableValue & _alpha_old;
    const VariableValue & _B_old; 
    const VariableValue & _xi_old; 
    const VariableValue & _I2_old;
    const VariableValue & _mu_old;
    const VariableValue & _lambda_old;
    const VariableValue & _gamma_old;
    /// mechanical strain rate
    const VariableValue & _mechanical_strain_rate;

    /// add options
    int _option;
    Real _Cd_constant;

    /// add healing option
    bool _healing;

    Real computeAlphaCr(Real xi);
};