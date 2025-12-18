#include "PoroStabilization.h"

registerMooseObject("farmsApp", PoroStabilization);

InputParameters
PoroStabilization::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription(
      "PSPG stabilization for equal-order u-p Biot poromechanics. "
      "Adds tau * (div(sigma_total)) * grad(q) to pressure equation.");
  
  params.addRequiredCoupledVar("displacements", "The displacement variables (disp_x, disp_y, disp_z)");
  
  return params;
}

PoroStabilization::PoroStabilization(const InputParameters & parameters)
  : Kernel(parameters),
    _ndisp(coupledComponents("displacements")),
    _disp_var(_ndisp),
    _tau_pspg(getMaterialProperty<Real>("tau_pspg")),
    _stress(getMaterialProperty<RankTwoTensor>("stress")),
    _Jacobian_mult(getMaterialProperty<RankFourTensor>("Jacobian_mult")),  // YOUR tangent!
    _stress_off_diag_jacobian(getMaterialProperty<RankTwoTensor>("stress_off_diag_jacobian"))  // YOUR ∂σ/∂p!
{
  // Couple displacement variables
  for (unsigned int i = 0; i < _ndisp; ++i)
  {
    _disp_var[i] = coupled("displacements", i);
  }
}

Real
PoroStabilization::computeQpResidual()
{
  // PSPG stabilization: tau_p * (div(sigma_total)) * grad(test_function)
  //
  // In weak form, div(σ) appears as: ∫ σ_ij * ∂ψ_i/∂x_j dΩ
  // We compute this contribution for stabilization
  
  RealVectorValue div_stress_vec;
  div_stress_vec.zero();
  
  for (unsigned int i = 0; i < _ndisp; ++i)
  {
    for (unsigned int j = 0; j < _ndisp; ++j)
    {
      // Contribution from stress component σ_ij in direction i
      div_stress_vec(i) += _stress[_qp](i, j) * _grad_test[_i][_qp](j);
    }
  }
  
  // PSPG term: τ * (div(σ) · ∇q)
  Real residual = _tau_pspg[_qp] * (div_stress_vec * _grad_test[_i][_qp]);
  
  return residual;
}

Real
PoroStabilization::computeQpJacobian()
{
  // d/dp of (tau * div(sigma_total) * grad(q))
  //
  // ∂(div(σ_total))/∂p is given by YOUR material's stress_off_diag_jacobian
  // This is: ∂σ_ij/∂p = -α_eff * δ_ij  (from your material)
  
  RealVectorValue div_dsigma_dp;
  div_dsigma_dp.zero();
  
  for (unsigned int i = 0; i < _ndisp; ++i)
  {
    for (unsigned int j = 0; j < _ndisp; ++j)
    {
      // ∂σ_ij/∂p from your material
      div_dsigma_dp(i) += _stress_off_diag_jacobian[_qp](i, j) * _grad_test[_i][_qp](j);
    }
  }
  
  // Chain rule: τ * (∂(div(σ))/∂p · ∇q) * ∂p/∂p_j
  return _tau_pspg[_qp] * (div_dsigma_dp * _grad_test[_i][_qp]) * _grad_phi[_j][_qp].norm();
}

Real
PoroStabilization::computeQpOffDiagJacobian(unsigned int jvar)
{
  // d/du_k of (tau * div(sigma_total) * grad(q))
  //
  // ∂σ_ij/∂u_k is given by YOUR material's Jacobian_mult tangent modulus
  // This is: ∂σ_ij/∂ε_kl * ∂ε_kl/∂u_k
  
  for (unsigned int comp = 0; comp < _ndisp; ++comp)
  {
    if (jvar == _disp_var[comp])
    {
      RealVectorValue div_dsigma_du;
      div_dsigma_du.zero();
      
      // For each stress component σ_ij
      for (unsigned int i = 0; i < _ndisp; ++i)
      {
        for (unsigned int j = 0; j < _ndisp; ++j)
        {
          // ∂σ_ij/∂u_comp via chain rule:
          // ∂σ_ij/∂u_comp = ∂σ_ij/∂ε_kl * ∂ε_kl/∂u_comp
          
          Real dsigma_ij_du_comp = 0.0;
          
          for (unsigned int k = 0; k < _ndisp; ++k)
          {
            for (unsigned int l = 0; l < _ndisp; ++l)
            {
              // Strain derivative: ∂ε_kl/∂u_comp
              // ε_kl = 0.5 * (∂u_k/∂x_l + ∂u_l/∂x_k)
              // ∂ε_kl/∂u_comp = 0.5 * (∂φ_j/∂x_l * δ_kc + ∂φ_j/∂x_k * δ_lc)
              
              Real deps_kl_du_comp = 0.0;
              if (k == comp)
                deps_kl_du_comp += 0.5 * _grad_phi[_j][_qp](l);
              if (l == comp)
                deps_kl_du_comp += 0.5 * _grad_phi[_j][_qp](k);
              
              // Use YOUR tangent modulus from material
              dsigma_ij_du_comp += _Jacobian_mult[_qp](i, j, k, l) * deps_kl_du_comp;
            }
          }
          
          // Accumulate divergence contribution
          div_dsigma_du(i) += dsigma_ij_du_comp * _grad_test[_i][_qp](j);
        }
      }
      
      // PSPG contribution: τ * (∂(div(σ))/∂u · ∇q)
      return _tau_pspg[_qp] * (div_dsigma_du * _grad_test[_i][_qp]);
    }
  }
  
  return 0.0;
}