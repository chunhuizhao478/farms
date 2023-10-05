/*
Update Breakage Variable Using Runge-Kutta own solver

- 10/5/2023 Chunhui Zhao

Include power-law correction on Cd (function of strain rate e)
if (e < 1e-4){ Cd = 10 };
else{ Cd = pow(10, log10(1+m*log10(e/1e-4)) ) }

Cb = CdCb_multiplier * Cd

*/

#pragma once

#include "AuxKernel.h"

class BreakageVarUpdate : public AuxKernel
{
    public:

    static InputParameters validParams();
    BreakageVarUpdate(const InputParameters & parameters);

    protected:

    virtual Real computeValue() override;

    /// constant parameters
    Real _Cd_min; //minimum Cd value for small strain (e < 1e-4)
    Real _C_BH;
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

    Real OdeIntegrator(Real alpha,
                       Real B,
                       Real I2,
                       Real xi,
                       Real mu,
                       Real gamma_damaged,
                       Real lambda,
                       Real C_B);
    Real computeBreakageEvolution(Real alpha,
                                  Real B,
                                  Real I2,
                                  Real xi,
                                  Real mu,
                                  Real gamma_damaged,
                                  Real lambda,
                                  Real C_B);
    Real computeGranularStateProb(Real alpha, Real xi);
    Real computeAlphaCr(Real xi);

};