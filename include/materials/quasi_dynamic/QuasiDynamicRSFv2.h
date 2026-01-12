#pragma once

#include "CZMComputeLocalTractionTotalBaseQDRSF2D.h"

class QuasiDynamicRSFv2 : public CZMComputeLocalTractionTotalBaseQDRSF2D
{
public:
    static InputParameters validParams();
    QuasiDynamicRSFv2(const InputParameters & parameters);

protected:
    void computeInterfaceTractionAndDerivatives() override;

private:
    // Rate-and-state friction parameters
    const Real _V_o;         // Reference slip rate
    const Real _f_o;         // Initial friction coefficient
    const Real _a;           // Direct effect parameter 
    const Real _b;           // State variable evolution parameter
    const Real _L;           // State variable characteristic distance

    // Radiation damping parameters
    const Real _shear_modulus;
    const Real _shear_wave_velocity;
    const Real _nu;          // Radiation damping coefficient

    // Enhanced weakening parameters
    const bool _enhanced_weakening;
    const Real _f_w;         // Weakened friction coefficient
    const Real _V_w;         // Weakening slip rate

    // Input variables
    const VariableValue& _normal_traction;
    const VariableValue& _tangential_traction;
    const VariableValue& _p_p;  // Neighbor pore pressure
    const VariableValue& _p_n;  // Local pore pressure

    // Background stress
    const Real _background_normal_stress;
    const Real _background_tangential_stress;

    // Old material properties
    const MaterialProperty<Real>& _sliprate_old;
    const MaterialProperty<Real>& _statevar_old;
    const MaterialProperty<Real>& _statevar_dot_old;
    const MaterialProperty<Real>& _slip_old;

    // Quasi-dynamic solver methods
    Real solveSlipRate(Real Tx, Real Ty, Real p_avg, Real theta_old, Real V_old);
    Real computeSlip(Real slip_old, Real V_new, Real V_old);
    Real computeStateVariable(Real theta_old, Real theta_dot, Real theta_dot_old);
};