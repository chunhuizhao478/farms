#pragma once
#include "InterfaceMaterial.h"
#include "ElasticDGUtils3D.h"

/*====================================================================
  FaultSlipWeakeningDG3DMaterial
  ------------------------------------------------------------------
  Centralizes fault interface physics for the slip‑weakening law:
    1. Rotate +/- side velocity & stress into local (n,t1,t2) frame.
    2. Form elastic Godunov (Riemann) solution for (un, σ_nn, u_t, τ_t).
    3. Apply linear slip‑weakening Coulomb criterion using STATEFUL
       accumulated slip magnitude (no external Aux slip variable).
    4. Compute imposed frictional traction τ~ (stick or sliding) and
       radiation damping adjusted tangential velocities u_t^{*}.
    5. Assemble numerical flux components in global coordinates.
    6. Update slip magnitude: d^{n+1} = d^{n} + |Δd|, with |Δd| = |slip_rate| * dt.

  Exposed Material Properties (all at qp on the interface):
    fault_flux_{sxx,sxy,sxz,syy,syz,szz,ux,uy,uz}  --> Consumed by DG kernel
    fault_tau_{t1,t2,mag}                          --> Frictional shear tractions
    fault_sigma_n                                  --> Effective normal stress (σ_n^G + σ0n)
    fault_mu_f                                     --> Current weakened friction coefficient μ(d)
    fault_slip_rate_{t1,t2,mag}                    --> Slip rate components / magnitude
    fault_ut{1,2}_star                             --> Star tangential velocities (for diagnostics)
    fault_slip                                     --> Accumulated slip magnitude (stateful)

  Notes:
    - Uses getMaterialPropertyOld to access previous slip for weakening; avoids
      need for an AuxVariable (user request).
    - Jacobian terms are currently handled in the DG kernel (analytical TBD).
====================================================================*/
class FaultSlipWeakeningDG3DMaterial : public InterfaceMaterial
{
public:
  static InputParameters validParams();
  FaultSlipWeakeningDG3DMaterial(const InputParameters & params);

protected:
  void initQpStatefulProperties() override; // zero slip at start
  void computeQpProperties() override;      // main physics

  // Coupled primal variables (minus / plus sides)
  const VariableValue & _ux; const VariableValue & _ux_n;
  const VariableValue & _uy; const VariableValue & _uy_n;
  const VariableValue & _uz; const VariableValue & _uz_n;
  const VariableValue & _sxx; const VariableValue & _sxx_n;
  const VariableValue & _syy; const VariableValue & _syy_n;
  const VariableValue & _szz; const VariableValue & _szz_n;
  const VariableValue & _sxy; const VariableValue & _sxy_n;
  const VariableValue & _sxz; const VariableValue & _sxz_n;
  const VariableValue & _syz; const VariableValue & _syz_n;

  // Elastic properties
  const MaterialProperty<Real> & _lambda;  const MaterialProperty<Real> & _lambda_n;
  const MaterialProperty<Real> & _mu;      const MaterialProperty<Real> & _mu_n;
  const MaterialProperty<Real> & _rho;     const MaterialProperty<Real> & _rho_n;

  // Friction & weakening constants
  const Real _mu_s_const; const Real _mu_d; const Real _Dc;
  const Real _tau0_t1_const; const Real _tau0_t2_const; const Real _sigma0n_const;
  // Optional spatial overrides (face values)
  const VariableValue * _mu_s_aux = nullptr;
  const VariableValue * _tau0_t1_aux = nullptr;
  const VariableValue * _sigma0n_aux = nullptr;

  // Stateful slip magnitude (current & old)
  MaterialProperty<Real> & _slip;                 // fault_slip (current step)
  const MaterialProperty<Real> & _slip_old;       // previous step

  // Flux properties (global frame)
  MaterialProperty<Real> & _flux_sxx; MaterialProperty<Real> & _flux_sxy; MaterialProperty<Real> & _flux_sxz;
  MaterialProperty<Real> & _flux_syy; MaterialProperty<Real> & _flux_syz; MaterialProperty<Real> & _flux_szz;
  MaterialProperty<Real> & _flux_ux;  MaterialProperty<Real> & _flux_uy;  MaterialProperty<Real> & _flux_uz;

  // Diagnostics
  MaterialProperty<Real> & _traction_t1; MaterialProperty<Real> & _traction_t2; MaterialProperty<Real> & _traction_mag;
  MaterialProperty<Real> & _normal_stress; MaterialProperty<Real> & _mu_f_prop;
  MaterialProperty<Real> & _sr_t1; MaterialProperty<Real> & _sr_t2; MaterialProperty<Real> & _sr_mag;
  MaterialProperty<Real> & _ut1_star; MaterialProperty<Real> & _ut2_star;

  inline Real weakeningMu(Real d, Real mu_s_val) const
  { return d < _Dc ? std::max(_mu_d, mu_s_val - (mu_s_val - _mu_d) * (d/_Dc)) : _mu_d; }
};
