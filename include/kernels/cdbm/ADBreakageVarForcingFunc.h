/*
Update Breakage Variable Using Runge-Kutta own solver

- 10/5/2023 Chunhui Zhao

Include power-law correction on Cd (function of strain rate e)
if (e < 1e-4){ Cd = 10 };
else{ Cd = pow(10, log10(1+m*log10(e/1e-4)) ) }

Cb = CdCb_multiplier * Cd

*/

#pragma once

#include "ADKernel.h"
#include "Material.h"

class ADBreakageVarForcingFunc : public ADKernel
{
public:
  static InputParameters validParams();

  ADBreakageVarForcingFunc(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual();

private:
    /// multiplier between Cb and Cbh
    ADReal _CBCBH_multiplier;
    ADReal _a0;
    ADReal _a1;
    ADReal _a2;
    ADReal _a3;
    ADReal _xi_min;
    ADReal _xi_max;
    ADReal _xi_d;
    ADReal _xi_0;
    ADReal _xi_1;
    ADReal _beta_width;
    /// multiplier between Cd and Cb
    ADReal _CdCb_multiplier;
    /// variable parameters
    const ADVariableValue & _alpha_old;
    const ADVariableValue & _B_old; 
    const ADVariableValue & _xi_old; 
    const ADVariableValue & _I2_old;
    const ADVariableValue & _mu_old;
    const ADVariableValue & _lambda_old;
    const ADVariableValue & _gamma_old;
    ADReal _Cd_constant;
    ADReal computeAlphaCr(ADReal xi);
};