//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov

#include "SmallStrainPlasticVolumetricStrainCoupling.h"

registerMooseObject("farmsApp", SmallStrainPlasticVolumetricStrainCoupling);

InputParameters
SmallStrainPlasticVolumetricStrainCoupling::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Volumetric coupling kernel for plastic deformation effects on pore pressure (small strain)");
  
  params.addRequiredCoupledVar("displacements", "The displacement variables");
  params.addParam<std::string>("base_name", "", "Optional parameter for multiple material systems");
  
  return params;
}

SmallStrainPlasticVolumetricStrainCoupling::SmallStrainPlasticVolumetricStrainCoupling(const InputParameters & parameters)
  : Kernel(parameters),
    _ndisp(coupledComponents("displacements")),
    _disp_var(_ndisp),
    _eps_p(getMaterialPropertyByName<RankTwoTensor>("eps_p")),
    _eps_p_old(getMaterialPropertyOldByName<RankTwoTensor>("eps_p")),
    _deps_p_dp(getMaterialProperty<RankTwoTensor>("deps_p_dp")),
    _deps_p_deps(getMaterialProperty<RankFourTensor>("deps_p_deps"))
{
  // Get displacement variable numbers for off-diagonal Jacobian
  for (unsigned int i = 0; i < _ndisp; ++i)
    _disp_var[i] = coupled("displacements", i);
}

Real
SmallStrainPlasticVolumetricStrainCoupling::computeQpResidual()
{
  // Volumetric plastic strain rate: tr(ε̇ᵖ) = [tr(εᵖ) - tr(εᵖ_old)] / dt
  Real tr_eps_p = _eps_p[_qp].trace();
  Real tr_eps_p_old = _eps_p_old[_qp].trace();
  Real tr_eps_p_dot = (tr_eps_p - tr_eps_p_old) / _dt;
  
  // Residual: ψ * tr(ε̇ᵖ)
  return _test[_i][_qp] * tr_eps_p_dot;
}

Real
SmallStrainPlasticVolumetricStrainCoupling::computeQpJacobian()
{
  // Jacobian: ∂R/∂p = ψ * φ * (1/dt) * ∂tr(εᵖ)/∂p
  //                  = ψ * φ * (1/dt) * tr(∂εᵖ/∂p)
  
  Real dtr_eps_p_dp = _deps_p_dp[_qp].trace();
  
  return _test[_i][_qp] * _phi[_j][_qp] * dtr_eps_p_dp / _dt;
}

Real
SmallStrainPlasticVolumetricStrainCoupling::computeQpOffDiagJacobian(unsigned int jvar)
{
  for (unsigned int comp = 0; comp < _ndisp; ++comp)
  {
    if (jvar == _disp_var[comp])
    {
      // ∂R/∂u = ψ * (1/dt) * ∂tr(εᵖ)/∂u
      // Chain rule: ∂tr(εᵖ)/∂u = Σ_k ∂εᵖ_kk/∂ε_mn * ∂ε_mn/∂u
      
      Real dtr_eps_p_du = 0.0;
      
      for (unsigned int k = 0; k < 3; ++k)          // trace over k
      {
        for (unsigned int m = 0; m < 3; ++m)        // strain indices
        {
          for (unsigned int n = 0; n < 3; ++n)
          {
            // Small strain: ∂ε_mn/∂u_comp = 1/2(δ_m,comp * ∂φ_j/∂x_n + δ_n,comp * ∂φ_j/∂x_m)
            Real deps_du = 0.5 * ((comp == m ? _grad_phi[_j][_qp](n) : 0.0) + 
                                  (comp == n ? _grad_phi[_j][_qp](m) : 0.0));
            
            dtr_eps_p_du += _deps_p_deps[_qp](k, k, m, n) * deps_du;
          }
        }
      }
      
      return _test[_i][_qp] * dtr_eps_p_du / _dt;
    }
  }
  
  return 0.0;
}