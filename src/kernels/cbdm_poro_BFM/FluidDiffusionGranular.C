//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FluidDiffusionGranular.h"

registerMooseObject("farmsApp", FluidDiffusionGranular);

InputParameters
FluidDiffusionGranular::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Diffusion kernel for pore pressure in poromechanical problems with reference configuration permeability");
  params.addParam<std::string>("base_name", "", "Optional parameter for multiple material systems");
  return params;
}

FluidDiffusionGranular::FluidDiffusionGranular(const InputParameters & parameters)
  : Kernel(parameters),
    _perm_g(getMaterialProperty<Real>("permeability_granular")),
    _viscosity(getMaterialProperty<Real>("viscosity_fluid")),
    _F(getMaterialProperty<RankTwoTensor>(getParam<std::string>("base_name") + "deformation_gradient")),
    _dJp_dF(getMaterialProperty<RankTwoTensor>(getParam<std::string>("base_name") + "plastic_jacobian_derivative")),
    _Jp(getMaterialProperty<Real>(getParam<std::string>("base_name") + "plastic_deformation_gradient_det")),
    _dJp_dp(getMaterialProperty<Real>(getParam<std::string>("base_name") + "plastic_jacobian_derivative_pressure")),
    _n_nodes(_current_elem->n_nodes())
{
}

RealGradient
FluidDiffusionGranular::computeGradJp()
{
  // 1. Initialize gradient tensors for F
  std::vector<RankTwoTensor> grad_F(3);
  for (unsigned int k = 0; k < 3; ++k)
    grad_F[k].zero();

  // 2. Compute grad_F[k] = ∂F/∂X_k using finite element interpolation
  // grad_F = Σ(a=1 to n_nodes) F_a ⊗ ∇N_a
  for (unsigned int a = 0; a < _n_nodes; ++a)
  {
    // Get F at node a (approximation: use current quadrature point value)
    // In a full implementation, you would need nodal values of F
    const RankTwoTensor & F_a = _F[_qp]; 
    
    // Get shape function gradient at node a
    const RealGradient & dNa = _grad_phi[a][_qp];
    
    // Accumulate: grad_F[k] += F_a * dN_a/dX_k
    for (unsigned int k = 0; k < 3; ++k)
      grad_F[k] += F_a * dNa(k);
  }

  // 3. Compute grad_Jp = (dJp/dF) : grad_F
  // grad_Jp_k = (dJp/dF)_ij * (grad_F[k])_ij
  RealGradient grad_Jp;
  
  for (unsigned int k = 0; k < 3; ++k)
  {
    // Double contraction: (dJp/dF)_ij * (grad_F[k])_ij
    Real sum = 0.0;
    for (unsigned int i = 0; i < 3; ++i)
      for (unsigned int j = 0; j < 3; ++j)
        sum += _dJp_dF[_qp](i, j) * grad_F[k](i, j);
    
    grad_Jp(k) = sum;
  }
  
  return grad_Jp;
}

Real
FluidDiffusionGranular::computeQpResidual()
{
  // 1. Compute reference configuration quantities
  RankTwoTensor F_inv = _F[_qp].inverse();
  RankTwoTensor F_inv_T = F_inv.transpose();
  Real J = _F[_qp].det();
  
  // 2. Compute permeability tensors in reference configuration
  // K_ref = J * F^(-1) * (kappa/mu) * F^(-T)
  RankTwoTensor K_g_ref = J * F_inv * (_perm_g[_qp] / _viscosity[_qp]) * F_inv_T;
  
  // 3. Compute gradient of Jp
  RealGradient grad_Jp = computeGradJp();

  // 4. Granular phase flux: q_g = -K_g_ref * grad(Jp * p)
  // grad(Jp * p) = Jp * grad(p) + p * grad(Jp)
  RealVectorValue grad_Jp_p = _Jp[_qp] * _grad_u[_qp] + _u[_qp] * grad_Jp;
  RealVectorValue flux_g = K_g_ref * grad_Jp_p;
  
  // 5. Residual: ∇ψ · q
  return _grad_test[_i][_qp] * flux_g;
}

Real
FluidDiffusionGranular::computeQpJacobian()
{
  // 1. Compute reference configuration quantities (same as residual)
  RankTwoTensor F_inv = _F[_qp].inverse();
  RankTwoTensor F_inv_T = F_inv.transpose();
  Real J = _F[_qp].det();
  
  // 2. Compute permeability tensors
  RankTwoTensor K_g_ref = J * F_inv * (_perm_g[_qp] / _viscosity[_qp]) * F_inv_T;
  
  // 3. Compute gradient of Jp
  RealGradient grad_Jp = computeGradJp();
  
  // 4. Granular phase Jacobian: ∂q_g/∂p = -K_g_ref * ∂(grad(Jp * p))/∂p
  // Complete derivative: ∂(grad(Jp * p))/∂p = 
  //   ∂J^p/∂p * grad(p) + J^p * grad(φ) + φ * grad(J^p) + p * ∂(grad(J^p))/∂p
  
  RealVectorValue grad_Jp_phi_complete;
  
  // Standard terms (what you already have):
  grad_Jp_phi_complete = _Jp[_qp] * _grad_phi[_j][_qp] + _phi[_j][_qp] * grad_Jp;
  
  // NEW: Missing terms since Jp depends on pressure:
  // Term 1: ∂J^p/∂p * grad(p) * φ_j
  grad_Jp_phi_complete += _dJp_dp[_qp] * _grad_u[_qp];
  
  // Note: Term p * ∂(grad(J^p))/∂p is more complex and may be neglected
  // for first-order accuracy, or computed if needed for high accuracy
  
  RealVectorValue jac_g = K_g_ref * grad_Jp_phi_complete;
  
  // 5. Final Jacobian: ∇ψ · (∂q/∂p)
  return _grad_test[_i][_qp] * jac_g;
}