/*
Update Breakage Variable Using Runge-Kutta own solver
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
    Real _C_B;
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
    /// variable parameters
    const VariableValue & _alpha_old;
    const VariableValue & _B_old; 
    const VariableValue & _xi_old; 
    const VariableValue & _I2_old;
    const VariableValue & _mu_old;
    const VariableValue & _lambda_old;
    const VariableValue & _gamma_old;

    Real OdeIntegrator(Real alpha,
                       Real B,
                       Real I2,
                       Real xi,
                       Real mu,
                       Real gamma_damaged,
                       Real lambda);
    Real computeBreakageEvolution(Real alpha,
                                  Real B,
                                  Real I2,
                                  Real xi,
                                  Real mu,
                                  Real gamma_damaged,
                                  Real lambda);
    Real computeGranularStateProb(Real alpha, Real xi);
    Real computeAlphaCr(Real xi);

};

