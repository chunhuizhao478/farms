#pragma once
#include "DGKernel.h"
#include "ElasticDGUtils3D.h"

/** 
 * Discontinuous Galerkin numerical flux for 3D elastic wave propagation using exact Godunov solver
 *
 * Implementation follows Pelties et al. (2012) "Three-dimensional dynamic rupture simulation 
 * with a high-order discontinuous Galerkin method on unstructured tetrahedral meshes"
 * Geophysical Journal International, doi:10.1111/j.1365-246X.2012.05626.x
 *
 * Key algorithmic components:
 * 1. Local coordinate transformation: Global (x,y,z) -> Interface (n,t1,t2) using orthonormal basis
 * 2. Material property usage: Each element uses its own λ, μ, ρ (no averaging across interface)
 * 3. Exact 1D Riemann solver: Three independent wave modes (P, S1, S2) per equation (13)
 * 4. Coordinate back-transformation: Interface fluxes rotated back to global components
 *
 * The Godunov flux provides upwind stabilization for hyperbolic elastic wave system:
 *   ∂σ/∂t + A_σ∇u = 0    (stress evolution with flux matrix A_σ)
 *   ∂u/∂t + A_u∇σ = 0    (velocity evolution with flux matrix A_u)
 *
 * Godunov flux formulation per equation (13):
 * u_n* = u_n^+ + (σ_nn^+ - σ_nn^-) / (Z_p^+ + Z_p^-)
 * σ_nn* = σ_nn^+ + Z_p^+ * (u_n^+ - u_n^-) * Z_p^- / (Z_p^+ + Z_p^-)
 * Similar formulation for tangential components with S-wave impedances Z_s
 *
 * Local frame flux formulation:
 * For stress: F_n(σ) = [(λ+2μ)u_n* on nn; λu_n* on t1t1,t2t2; μu_tk* on ntk]  
 * For velocity: F_n(u) = [σ_nn/ρ, σ_nt1/ρ, σ_nt2/ρ] (Godunov stresses over element density)
 */
class ElasticGodunovDGFlux3D : public DGKernel
{
public:
  static InputParameters validParams();
  ElasticGodunovDGFlux3D(const InputParameters & params);

protected:
  Real computeQpResidual(Moose::DGResidualType type) override;
  // Jacobian contributions (currently omitted: return 0 to satisfy abstract interface)
  Real computeQpJacobian(Moose::DGJacobianType type) override;
  Real computeQpOffDiagJacobian(Moose::DGJacobianType type, unsigned int jvar) override;

  // Coupled variables (element and neighbor)
  const VariableValue & _ux; const VariableValue & _ux_neighbor;
  const VariableValue & _uy; const VariableValue & _uy_neighbor;
  const VariableValue & _uz; const VariableValue & _uz_neighbor;
  const VariableValue & _sxx; const VariableValue & _sxx_neighbor;
  const VariableValue & _syy; const VariableValue & _syy_neighbor;
  const VariableValue & _szz; const VariableValue & _szz_neighbor;
  const VariableValue & _sxy; const VariableValue & _sxy_neighbor;
  const VariableValue & _sxz; const VariableValue & _sxz_neighbor;
  const VariableValue & _syz; const VariableValue & _syz_neighbor;

  // Material properties (both sides)
  const MaterialProperty<Real> & _lambda;  const MaterialProperty<Real> & _lambda_n;
  const MaterialProperty<Real> & _mu;      const MaterialProperty<Real> & _mu_n;
  const MaterialProperty<Real> & _rho;     const MaterialProperty<Real> & _rho_n;

  // Optional: skip flux on an interface between two specific subdomain IDs (fault)
  bool _use_skip_pair = false;
  SubdomainID _skip_a = Moose::INVALID_BLOCK_ID;
  SubdomainID _skip_b = Moose::INVALID_BLOCK_ID;
};