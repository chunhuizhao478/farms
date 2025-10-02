#pragma once
#include "DGKernel.h"

/** Fault DG flux with linear slip-weakening:
 *  1) Build Godunov state in (n,t1,t2).                                    [2,3]
 *  2) Test Coulomb failure |τ^G + τ0| ? μ(d) |σ_n^G + σ0n|.                 [2]
 *  3) If active, impose τ~_t = μ(d) |σ_n| * dir(τ^G + τ0); else τ~_t = τ^G. [2]
 *  4) Radiation damping: Δ\dot d_t = 2/Zs (τ~ - τ^G); side updates u_t^{±,~}=u_t^± ± (τ~ - τ^±)/Zs. [3]
 *  5) Build flux with u_n^G, u_t^*, σ_nn^G, τ~_t and rotate back.
 */
class FaultSlipWeakeningDG3D : public DGKernel {
public:
  static InputParameters validParams();
  FaultSlipWeakeningDG3D(const InputParameters & params);
protected:
  Real computeQpResidual(Moose::DGResidualType type) override;
  Real computeQpJacobian(Moose::DGJacobianType type) override;
  Real computeQpOffDiagJacobian(Moose::DGJacobianType type, unsigned int jvar) override;
  // Flux material properties supplied by FaultSlipWeakeningDG3DMaterial
  const MaterialProperty<Real> & _fx_sxx; const MaterialProperty<Real> & _fx_sxy; const MaterialProperty<Real> & _fx_sxz;
  const MaterialProperty<Real> & _fx_syy; const MaterialProperty<Real> & _fx_syz; const MaterialProperty<Real> & _fx_szz;
  const MaterialProperty<Real> & _fx_ux;  const MaterialProperty<Real> & _fx_uy;  const MaterialProperty<Real> & _fx_uz;
};