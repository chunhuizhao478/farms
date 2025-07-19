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
    _F(getMaterialProperty<RankTwoTensor>(getParam<std::string>("base_name") + "deformation_gradient")),
    _Jp(getMaterialProperty<Real>(getParam<std::string>("base_name") + "plastic_deformation_gradient_det")),
    _dJp_dp(getMaterialProperty<Real>(getParam<std::string>("base_name") + "plastic_jacobian_derivative_pressure"))
{
}

Real
FluidDiffusionGranular::computeQpResidual()
{
  // 1. Compute reference configuration quantities
  RankTwoTensor F_inv = _F[_qp].inverse();
  RankTwoTensor F_inv_T = F_inv.transpose();
  Real J = _F[_qp].det();
  
  // 2. Compute permeability tensor in reference configuration
  // Option A: Full reference configuration permeability
  // RankTwoTensor K_g_ref = J * F_inv * (_perm_g[_qp] / _viscosity[_qp]) * F_inv_T;
  RankTwoTensor K_g_ref = (_perm_g[_qp]) * RankTwoTensor::Identity();

  // 3. Compute gradient of Jp using pressure dependence approximation
  // ∇J^p ≈ (∂J^p/∂p) * ∇p
  RealGradient grad_Jp = _dJp_dp[_qp] * _grad_u[_qp];
  
  // 4. Granular phase flux: q_g = -K_g_ref * grad(J^p * p)
  // grad(J^p * p) = J^p * grad(p) + p * grad(J^p)
  RealVectorValue grad_Jp_p = _Jp[_qp] * _grad_u[_qp] + _u[_qp] * grad_Jp;
  RealVectorValue flux_g = K_g_ref * grad_Jp_p;
  
  // 5. Residual: ∇ψ · q_g
  return _grad_test[_i][_qp] * flux_g;
}

Real
FluidDiffusionGranular::computeQpJacobian()
{
  // 1. Compute reference configuration quantities (same as residual)
  RankTwoTensor F_inv = _F[_qp].inverse();
  RankTwoTensor F_inv_T = F_inv.transpose();
  Real J = _F[_qp].det();
  
  // 2. Compute permeability tensor in reference configuration
  // Option A: Full reference configuration permeability
  // RankTwoTensor K_g_ref = J * F_inv * (_perm_g[_qp] / _viscosity[_qp]) * F_inv_T;
  RankTwoTensor K_g_ref = (_perm_g[_qp]) * RankTwoTensor::Identity();

  // 3. Compute gradient of Jp (same as residual)
  RealGradient grad_Jp = _dJp_dp[_qp] * _grad_u[_qp];
  
  // 4. Jacobian: ∂q_g/∂p = -K_g_ref * ∂(grad(J^p * p))/∂p
  // 
  // Complete derivative: ∂(grad(J^p * p))/∂p = 
  //   ∂J^p/∂p * grad(p) * φ_j + J^p * grad(φ_j) + φ_j * grad(J^p) + p * ∂(grad(J^p))/∂p
  //
  // Where: ∂(grad(J^p))/∂p = ∂/∂p[∂J^p/∂p * grad(p)] = ∂J^p/∂p * grad(φ_j)
  
  RealVectorValue jac_grad_Jp_p;
  
  // Term 1: ∂J^p/∂p * grad(p) * φ_j
  jac_grad_Jp_p += _dJp_dp[_qp] * _grad_u[_qp] * _phi[_j][_qp];
  
  // Term 2: J^p * grad(φ_j)
  jac_grad_Jp_p += _Jp[_qp] * _grad_phi[_j][_qp];
  
  // Term 3: φ_j * grad(J^p)
  jac_grad_Jp_p += _phi[_j][_qp] * grad_Jp;
  
  // Term 4: p * ∂(grad(J^p))/∂p = p * ∂J^p/∂p * grad(φ_j)
  jac_grad_Jp_p += _u[_qp] * _dJp_dp[_qp] * _grad_phi[_j][_qp];
  
  RealVectorValue jac_flux = K_g_ref * jac_grad_Jp_p;
  
  // 5. Final Jacobian: ∇ψ · (∂q_g/∂p)
  return _grad_test[_i][_qp] * jac_flux;
}