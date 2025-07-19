//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PlasticVolumetricStrainCoupling.h"

registerMooseObject("farmsApp", PlasticVolumetricStrainCoupling);

InputParameters
PlasticVolumetricStrainCoupling::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Volumetric coupling kernel for plastic deformation effects on pore pressure");
  
  params.addRequiredCoupledVar("displacements", "The displacement variables");
  params.addParam<std::string>("base_name", "", "Optional parameter for multiple material systems");
  
  return params;
}

PlasticVolumetricStrainCoupling::PlasticVolumetricStrainCoupling(const InputParameters & parameters)
  : Kernel(parameters),
    _ndisp(coupledComponents("displacements")),
    _disp_var(_ndisp),
    _Jp(getMaterialProperty<Real>(getParam<std::string>("base_name") + "plastic_deformation_gradient_det")),
    _Dp(getMaterialProperty<RankTwoTensor>(getParam<std::string>("base_name") + "plastic_strain_rate")),
    _dJp_dF(getMaterialProperty<RankTwoTensor>(getParam<std::string>("base_name") + "plastic_jacobian_derivative")),
    _dJp_dp(getMaterialProperty<Real>(getParam<std::string>("base_name") + "plastic_jacobian_derivative_pressure")),
    _dDp_dF(getMaterialProperty<RankFourTensor>(getParam<std::string>("base_name") + "plastic_deformation_rate_derivative")),
    _dDp_dp(getMaterialProperty<RankTwoTensor>(getParam<std::string>("base_name") + "plastic_deformation_rate_derivative_pressure"))
{
  // Get displacement variable numbers for off-diagonal Jacobian
  for (unsigned int i = 0; i < _ndisp; ++i)
    _disp_var[i] = coupled("displacements", i);
}

Real
PlasticVolumetricStrainCoupling::computeQpResidual()
{
  // Compute trace of plastic strain rate: tr(D^p)
  Real tr_Dp = _Dp[_qp].trace();
  
  // Residual: ψ * J^p * tr(D^p)
  return _test[_i][_qp] * _Jp[_qp] * tr_Dp;
}

Real
PlasticVolumetricStrainCoupling::computeQpJacobian()
{
  // Now we account for the dependence of J^p and D^p on pore pressure
  // Residual: R = ψ * J^p * tr(D^p)
  // Jacobian: ∂R/∂p = ψ * [∂J^p/∂p * tr(D^p) + J^p * ∂tr(D^p)/∂p]
  
  // Get current values
  Real tr_Dp = _Dp[_qp].trace();
  
  // Term 1: ∂J^p/∂p * tr(D^p)
  Real dtr_Dp_dp = _dDp_dp[_qp].trace();  
  Real dJp_dp = _dJp_dp[_qp];     
  
  // Term 2: J^p * ∂tr(D^p)/∂p
  // This is just: J^p * tr(∂D^p/∂p) = J^p * dtr_Dp_dp
  
  // Final Jacobian: ∂R/∂p = ψ * [∂J^p/∂p * tr(D^p) + J^p * ∂tr(D^p)/∂p]
  return _test[_i][_qp] * _phi[_j][_qp] * (dJp_dp * tr_Dp + _Jp[_qp] * dtr_Dp_dp);
}

Real
PlasticVolumetricStrainCoupling::computeQpOffDiagJacobian(unsigned int jvar)
{
  for (unsigned int i = 0; i < _ndisp; ++i)
  {
    if (jvar == _disp_var[i])
    {
      const Real tr_Dp = _Dp[_qp].trace();

      // dJp/du: chain rule ∂Jp/∂F_kl * ∂F_kl/∂u_j
      Real dJp_du = 0.0;
      for (unsigned int k = 0; k < 3; ++k)
        for (unsigned int l = 0; l < 3; ++l)
        {
          const Real dF_kl_du = (i == k) ? _grad_phi[_j][_qp](l) : 0.0;
          dJp_du += _dJp_dF[_qp](k, l) * dF_kl_du;
        }

      // d(tr(Dp))/du = trace of dDp/du = sum_k dDp_kk/du
      Real dtr_Dp_du = 0.0;
      for (unsigned int k = 0; k < 3; ++k)
        for (unsigned int m = 0; m < 3; ++m)
          for (unsigned int n = 0; n < 3; ++n)
          {
            const Real dF_mn_du = (i == m) ? _grad_phi[_j][_qp](n) : 0.0;
            dtr_Dp_du += _dDp_dF[_qp](k, k, m, n) * dF_mn_du;
          }

      // Final ∂R/∂u_j = ψ * [∂Jp/∂u_j * tr(Dp) + Jp * ∂tr(Dp)/∂u_j]
      return _test[_i][_qp] * (dJp_du * tr_Dp + _Jp[_qp] * dtr_Dp_du);
    }
  }

  return 0.0;
}