//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FluidDiffusion2.h"

registerMooseObject("farmsApp", FluidDiffusion2);

InputParameters
FluidDiffusion2::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Diffusion kernel for pore pressure in poromechanical problems with reference configuration permeability");
  params.addParam<std::string>("base_name", "", "Optional parameter for multiple material systems");
  params.addParam<bool>(
    "large_kinematics", false, "Set to true to use large kinematics (F, J, etc.) in the formulation");
  return params;
}

FluidDiffusion2::FluidDiffusion2(const InputParameters & parameters)
  : Kernel(parameters),
    _perm_s(getMaterialProperty<Real>("permeability_solid")),
    _F(getMaterialProperty<RankTwoTensor>(getParam<std::string>("base_name") + "deformation_gradient")),
    _large_kinematics(getParam<bool>("large_kinematics")) 
{
}

Real
FluidDiffusion2::computeQpResidual()
{
  // 1. Compute reference configuration quantities
  RankTwoTensor F_inv = _F[_qp].inverse();
  RankTwoTensor F_inv_T = F_inv.transpose();
  Real J = _F[_qp].det();
  
  // // 2. Compute permeability tensors in reference configuration
  // // K_ref = J * F^(-1) * (kappa/mu) * F^(-T)
  RankTwoTensor K_s_ref;
  // if (_large_kinematics){
  //    K_s_ref = J * F_inv * _perm_s[_qp]  * F_inv_T;
  // }
  // else{
  //    K_s_ref = _perm_s[_qp] * RankTwoTensor::Identity();
  // }

  K_s_ref = _perm_s[_qp] * RankTwoTensor::Identity();
 
  // 3. Compute flux components
  // Solid phase flux: q_s = K_s_ref * grad(p)
  RealVectorValue flux_s = K_s_ref * _grad_u[_qp];
  
  // 4. Residual: ∇ψ · q
  return _grad_test[_i][_qp] * flux_s;
}

Real
FluidDiffusion2::computeQpJacobian()
{
  // 1. Compute reference configuration quantities (same as residual)
  RankTwoTensor F_inv = _F[_qp].inverse();
  RankTwoTensor F_inv_T = F_inv.transpose();
  Real J = _F[_qp].det();
  
  // // 2. Compute permeability tensors
  // RankTwoTensor K_s_ref;
  // if (_large_kinematics){
  //    K_s_ref = J * F_inv * _perm_s[_qp]  * F_inv_T;
  // }
  // else{
  //    K_s_ref = _perm_s[_qp] * RankTwoTensor::Identity();
  // }
  RankTwoTensor K_s_ref;
  K_s_ref = _perm_s[_qp] * RankTwoTensor::Identity();

  // 3. Compute Jacobian components
  // Solid phase Jacobian: ∂q_s/∂p = K_s_ref * grad(φ)
  RealVectorValue jac_s = K_s_ref * _grad_phi[_j][_qp];
  
  // 4. Final Jacobian: ∇ψ · (∂q/∂p)
  return _grad_test[_i][_qp] * jac_s;
}