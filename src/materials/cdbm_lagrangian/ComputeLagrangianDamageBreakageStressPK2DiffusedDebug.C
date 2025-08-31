//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ComputeLagrangianDamageBreakageStressPK2DiffusedDebug.h"

registerMooseObject("farmsApp", ComputeLagrangianDamageBreakageStressPK2DiffusedDebug);

InputParameters
ComputeLagrangianDamageBreakageStressPK2DiffusedDebug::validParams()
{
  InputParameters params = ComputeLagrangianStressPK1::validParams();

  // --- NEW: controls for increment-based flip handling ---
  params.addParam<bool>("enforce_monotonic_no_negative_elastic", true,
                        "Under monotonic loading (based on total deviatoric strain increment), "
                        "project elastic deviatoric strain to zero instead of letting it flip sign; "
                        "assign the deviatoric total strain to plastic for that step.");
  params.addParam<Real>("flip_tol", 1e-12,
                        "Tolerance used in the deviatoric increment sign tests.");

  return params;
}

ComputeLagrangianDamageBreakageStressPK2DiffusedDebug::ComputeLagrangianDamageBreakageStressPK2DiffusedDebug(const InputParameters & parameters)
  : ComputeLagrangianStressPK1(parameters),
    _Fp(declareProperty<RankTwoTensor>(_base_name + "plastic_deformation_gradient")),
    _Jp(declareProperty<Real>(_base_name + "plastic_deformation_gradient_det")),
    _Fe(declareProperty<RankTwoTensor>(_base_name + "elastic_deformation_gradient")),
    _Tau(declareProperty<RankTwoTensor>(_base_name + "deviatroic_stress")),
    _Ee(declareProperty<RankTwoTensor>(_base_name + "green_lagrange_elastic_strain")),
    _Ep(declareProperty<RankTwoTensor>(_base_name + "plastic_strain")),
    _E(declareProperty<RankTwoTensor>(_base_name + "total_lagrange_strain")),
    _I1(declareProperty<Real>(_base_name + "first_elastic_strain_invariant")),
    _I2(declareProperty<Real>(_base_name + "second_elastic_strain_invariant")),
    _xi(declareProperty<Real>(_base_name + "strain_invariant_ratio")),
    _S(declareProperty<RankTwoTensor>(_base_name + "pk2_stress")),
    _Tp(declareProperty<RankTwoTensor>(_base_name + "plastic_stress")),
    _C(declareProperty<RankFourTensor>(_base_name + "pk2_jacobian")),
    _Dp(declareProperty<RankTwoTensor>(_base_name + "plastic_strain_rate")),
    _Fp_dot(declareProperty<RankTwoTensor>(_base_name + "cdbm_plastic_deformation_gradient_rate")),
    _F_dot(declareProperty<RankTwoTensor>(_base_name + "cdbm_deformation_gradient_rate")),
    _D(declareProperty<RankTwoTensor>(_base_name + "deformation_rate")),
    //---------------------------------------------------------------------------------------------//
    // Add option to add dilatancy/compaction effect //Follow paper Section 7.1
    _eta(declareProperty<Real>(_base_name + "plastic_volume_change")),
    _dilatancy_function_beta(declareProperty<Real>(_base_name + "dilatancy_function_beta")),
    _shear_rate_nu(declareProperty<Real>(_base_name + "shear_rate_nu")),
    //---------------------------------------------------------------------------------------------//
    _deviatroic_strain_rate(declareProperty<Real>("deviatroic_strain_rate")),
    //---------------------------------------------------------------------------------------------//
    _lambda_const(getMaterialProperty<Real>("lambda_const")),
    _shear_modulus(getMaterialProperty<Real>("shear_modulus")),
    _damaged_modulus(getMaterialProperty<Real>("damaged_modulus")),
    _B_breakagevar(getMaterialProperty<Real>("B_damagedvar")),
    _B_breakagevar_old(getMaterialPropertyOldByName<Real>("B_damagedvar")),
    _Tau_old(getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "deviatroic_stress")),
    _Fp_old(getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "plastic_deformation_gradient")),
    _F_old(getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "deformation_gradient")),
    _Ep_old(getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "plastic_strain")),
    //---------------------------------------------------------------------------------------------//
    // Add option to add dilatancy/compaction effect //Follow paper Section 7.1
    _eta_old(getMaterialPropertyOldByName<Real>(_base_name + "plastic_volume_change")),
    //---------------------------------------------------------------------------------------------//
    _C_g(getMaterialProperty<Real>("C_g")),
    _m1(getMaterialProperty<Real>("m1")),
    _m2(getMaterialProperty<Real>("m2")),
    _a0(getMaterialProperty<Real>("a0")),
    _a1(getMaterialProperty<Real>("a1")),
    _a2(getMaterialProperty<Real>("a2")),
    _a3(getMaterialProperty<Real>("a3")),
    //---------------------------------------------------------------------------------------------//
    _structural_stress_coefficient(getMaterialProperty<Real>("structural_stress_coefficient")),
    _grad_alpha_damagedvar(getMaterialProperty<RealGradient>("gradient_alpha_damagedvar")),
    //---------------------------------------------------------------------------------------------//
    //For the strain rate Cd
    _velgrad_L(getMaterialProperty<RankTwoTensor>("velgrad_L")),
    //---------------------------------------------------------------------------------------------//
    //use state dependent variables
    _Theta(declareProperty<Real>("state_variable")),
    _Theta_old(getMaterialPropertyOldByName<Real>("state_variable")),
    _use_state_var_evolution_mat(getMaterialProperty<bool>("use_state_var_evolution_mat")),
    _const_A_mat(getMaterialProperty<Real>("const_A_mat")),
    _const_B_mat(getMaterialProperty<Real>("const_B_mat")),
    _const_theta_o_mat(getMaterialProperty<Real>("const_theta_o_mat")),
    _initial_theta0_mat(getMaterialProperty<Real>("initial_theta0_mat")),
    //---------------------------------------------------------------------------------------------//
    //add shear stress perturbation
    _shear_stress_perturbation(getMaterialPropertyOldByName<Real>("shear_stress_perturbation"))
{
}

//Rules:See https://github.com/idaholab/moose/discussions/19450
//Only the object that declares the material property can assign values to it.
//Objects can request material properties, gaining read-only access to their values.
//When any object (including the object that declares it) requests the old value of a material property, that property becomes "stateful".
//All stateful material properties must be initialized within the initQpStatefulProperties call.
//
void
ComputeLagrangianDamageBreakageStressPK2DiffusedDebug::initQpStatefulProperties()
{
  _Fp[_qp] = RankTwoTensor::Identity();
  _Fe[_qp] = RankTwoTensor::Identity();
  _Tau[_qp].zero();
  _Ee[_qp].zero();
  _Ep[_qp].zero();
  _E[_qp].zero();
  _I1[_qp] = 0.0;
  _I2[_qp] = 0.0;
  _xi[_qp] = -sqrt(3);
  _S[_qp].zero();
  _C[_qp].zero();
  _eta[_qp] = 0.0;
  _dilatancy_function_beta[_qp] = 0.0;
  _shear_rate_nu[_qp] = 0.0;
  _Theta[_qp] = _initial_theta0_mat[_qp];

  _Fp_dot[_qp].zero();
  _F_dot[_qp].zero();
}

void
ComputeLagrangianDamageBreakageStressPK2DiffusedDebug::computeQpPK1Stress()
{
  //--------------------------------------------------------------------------
  // PK2 update (fills _S, _C, _Fp, _Fe, _Ee, etc.)
  computeQpPK2Stress();

  _Jp[_qp] = _Fp[_qp].det();

  RankTwoTensor Fpinv = _Fp[_qp].inverse();

  //--------------------------------------------------------------------------
  // Compute Fp_dot, F_dot (first-order rate using D and Dp)
  RankTwoTensor Fp_dot; Fp_dot.zero();
  RankTwoTensor F_dot;  F_dot.zero();

  for (unsigned int i = 0; i < 3; i++)
    for (unsigned int j = 0; j < 3; j++)
      for (unsigned int m = 0; m < 3; m++)
      {
        F_dot(i,j)  += _D[_qp](i,m)  * _F[_qp](m,j);
        Fp_dot(i,j) += _Dp[_qp](i,m) * _Fp[_qp](m,j);
      }

  // --- NEW: store into stateful properties used below
  _Fp_dot[_qp] = Fp_dot;
  _F_dot[_qp]  = F_dot;

  //--------------------------------------------------------------------------
  // Precompute dFpdF tensor
  Real dFpdF_tensor[3][3][3][3];
  for (unsigned int i = 0; i < 3; i++)
    for (unsigned int j = 0; j < 3; j++)
      for (unsigned int k = 0; k < 3; k++)
        for (unsigned int l = 0; l < 3; l++)
          dFpdF_tensor[i][j][k][l] =
              (_dt == 0.0 || _Fp_dot[_qp](i,j) == 0.0 || _F_dot[_qp](k,l) == 0.0)
              ? 0.0 : _Fp_dot[_qp](i,j) / _F_dot[_qp](k,l);

  auto delta = [](int i, int j) -> Real { return (i == j) ? 1.0 : 0.0; };

  //--------------------------------------------------------------------------
  // dJpdF
  Real dJpdF_tensor[3][3];
  for (unsigned int k = 0; k < 3; k++)
    for (unsigned int l = 0; l < 3; l++)
    {
      Real sum = 0.0;
      for (unsigned int m = 0; m < 3; m++)
        for (unsigned int n = 0; n < 3; n++)
          sum += _Jp[_qp] * Fpinv(n, m) * dFpdF_tensor[m][n][k][l];
      dJpdF_tensor[k][l] = sum;
    }

  //--------------------------------------------------------------------------
  // dFedF
  Real dFedF_tensor[3][3][3][3];
  for (unsigned int i = 0; i < 3; i++)
    for (unsigned int m = 0; m < 3; m++)
      for (unsigned int k = 0; k < 3; k++)
        for (unsigned int l = 0; l < 3; l++)
        {
          Real val = delta(i, k) * Fpinv(l, m);
          for (unsigned int h = 0; h < 3; h++)
            for (unsigned int r = 0; r < 3; r++)
              val -= _Fe[_qp](i, h) * dFpdF_tensor[h][r][k][l] * Fpinv(r, m);
          dFedF_tensor[i][m][k][l] = val;
        }

  //--------------------------------------------------------------------------
  // dEdF
  Real dEdF_tensor[3][3][3][3];
  for (unsigned int p = 0; p < 3; p++)
    for (unsigned int q = 0; q < 3; q++)
      for (unsigned int k = 0; k < 3; k++)
        for (unsigned int l = 0; l < 3; l++)
        {
          Real sum = 0.0;
          for (unsigned int m = 0; m < 3; m++)
            sum += 0.5 * ( dFedF_tensor[m][p][k][l] * _Fe[_qp](m, q)
                         + _Fe[_qp](m, p) * dFedF_tensor[m][q][k][l] );
          dEdF_tensor[p][q][k][l] = sum;
        }

  //--------------------------------------------------------------------------
  // dFpmdF
  Real dFpmdF_tensor[3][3][3][3];
  for (unsigned int j = 0; j < 3; j++)
    for (unsigned int n = 0; n < 3; n++)
      for (unsigned int k = 0; k < 3; k++)
        for (unsigned int l = 0; l < 3; l++)
        {
          Real sum = 0.0;
          for (unsigned int i = 0; i < 3; i++)
            for (unsigned int m = 0; m < 3; m++)
              sum += - Fpinv(j, i) * dFpdF_tensor[i][m][k][l] * Fpinv(m, n);
          dFpmdF_tensor[j][n][k][l] = sum;
        }

  //--------------------------------------------------------------------------
  // PK1 jacobian
  RankFourTensor pk_jacobian_val;
  pk_jacobian_val.zero();

  for (unsigned int i = 0; i < 3; i++)
    for (unsigned int j = 0; j < 3; j++)
      for (unsigned int k = 0; k < 3; k++)
        for (unsigned int l = 0; l < 3; l++)
        {
          Real accum = 0.0;
          // Term 1
          for (unsigned int m = 0; m < 3; m++)
            for (unsigned int n = 0; n < 3; n++)
              accum += dJpdF_tensor[k][l] * _Fe[_qp](i, m) * _S[_qp](m, n) * Fpinv(j, n);
          // Term 2
          for (unsigned int m = 0; m < 3; m++)
            for (unsigned int n = 0; n < 3; n++)
              accum += _Jp[_qp] * dFedF_tensor[i][m][k][l] * _S[_qp](m, n) * Fpinv(j, n);
          // Term 3
          for (unsigned int m = 0; m < 3; m++)
            for (unsigned int n = 0; n < 3; n++)
              for (unsigned int p = 0; p < 3; p++)
                for (unsigned int q = 0; q < 3; q++)
                  accum += _Jp[_qp] * _Fe[_qp](i, m) * _C[_qp](m, n, p, q) * dEdF_tensor[p][q][k][l] * Fpinv(j, n);
          // Term 4
          for (unsigned int m = 0; m < 3; m++)
            for (unsigned int n = 0; n < 3; n++)
              accum += _Jp[_qp] * _Fe[_qp](i, m) * _S[_qp](m, n) * dFpmdF_tensor[j][n][k][l];

          pk_jacobian_val(i, j, k, l) = accum;
        }

  //--------------------------------------------------------------------------
  // PK2-to-PK1 wrapping: using large kinematics formulation.
  if (_large_kinematics)
  {
    _pk1_stress[_qp] = _Jp[_qp] * _Fe[_qp] * _S[_qp] * Fpinv.transpose();
    _pk1_jacobian[_qp] = pk_jacobian_val;
  }
  else
    mooseError("Must select 'large_kinematics' option!");
}

void
ComputeLagrangianDamageBreakageStressPK2DiffusedDebug::computeQpPK2Stress()
{
  // --- plastic update ---
  RankTwoTensor Fp_updated = computeQpFp();

  const RankTwoTensor I = RankTwoTensor::Identity();

  // current kinematics
  RankTwoTensor Fe = _F[_qp] * Fp_updated.inverse();
  RankTwoTensor Ee_trial = 0.5 * (Fe.transpose() * Fe - I);                 // trial elastic
  RankTwoTensor Etot_new = 0.5 * (_F[_qp].transpose() * _F[_qp] - I);       // total GL strain

  // old kinematics
  RankTwoTensor Fe_old = _F_old[_qp] * _Fp_old[_qp].inverse();
  RankTwoTensor Ee_old = 0.5 * (Fe_old.transpose() * Fe_old - I);
  RankTwoTensor Etot_old = 0.5 * (_F_old[_qp].transpose() * _F_old[_qp] - I);

  auto deviator = [&](const RankTwoTensor & A) { return A - (A.trace()/3.0) * I; };
  auto dot4 = [](const RankTwoTensor & A, const RankTwoTensor & B){
    Real s = 0.0;
    for (unsigned i=0;i<3;++i)
      for (unsigned j=0;j<3;++j)
        s += A(i,j) * B(i,j);
    return s;
  };

  // increment-based test (deviatoric)
  RankTwoTensor dEtot_dev   = deviator(Etot_new - Etot_old);
  RankTwoTensor Ee_dev_old  = deviator(Ee_old);
  RankTwoTensor Ee_dev_tr   = deviator(Ee_trial);

  const Real flip_tol  = getParam<Real>("flip_tol");
  const bool enforce   = getParam<bool>("enforce_monotonic_no_negative_elastic");

  const Real s_monotonic = dot4(dEtot_dev, Ee_dev_old);   // >0 ⇒ monotonic shear increment
  const Real s_trial     = dot4(Ee_dev_tr, Ee_dev_old);   // ≤0 ⇒ elastic shear would flip

  const bool old_has_shear = (std::sqrt(dot4(Ee_dev_old, Ee_dev_old)) > 1e-16);
  const bool monotonic_inc = (s_monotonic >  flip_tol);
  const bool would_flip    = (s_trial     <= flip_tol);

  // initialize with trial
  RankTwoTensor Ee = Ee_trial;

  RankTwoTensor Ep_fromFp = 0.5 * (Fp_updated.transpose() * Fp_updated - I); // plastic (trial)
  RankTwoTensor Ep = Ep_fromFp;                                              // may be corrected
  RankTwoTensor E  = Fp_updated.transpose() * Ee * Fp_updated + Ep;          // may be corrected

  // If monotonic loading and elastic deviatoric would flip, enforce:
  //   Ee_dev -> 0  and  Ep_dev = dev(E_tot_new)
  if (enforce && old_has_shear && monotonic_inc && would_flip)
  {
    // project elastic to pure volumetric
    const Real I1_tr = Ee_trial.trace();
    Ee = (I1_tr/3.0) * I;

    // plastic deviatoric takes all total deviatoric
    RankTwoTensor Ep_dev_new = deviator(Etot_new);

    // keep plastic volumetric part from trial (or your dilatancy model)
    const Real Ep_vol = (Ep_fromFp.trace())/3.0;

    Ep = Ep_dev_new + Ep_vol * I;

    // rebuild total strain with corrected split
    E  = Fp_updated.transpose() * Ee * Fp_updated + Ep;
  }

  // invariants from Ee actually used
  Real I1 = Ee.trace();
  Real I2 = 0.0;
  for (unsigned int i = 0; i < 3; ++i)
    for (unsigned int j = 0; j < 3; ++j)
      I2 += Ee(i,j) * Ee(i,j);

  // xi (guard small I2)
  const Real adjusted_I2 = std::max(I2, 1e-12);
  const Real sqrt_I2     = std::sqrt(adjusted_I2);
  Real xi = I1 / sqrt_I2;
  if (std::isnan(xi))
    xi = -std::sqrt(3.0);

  // PK2 stress (solid/granular mix)
  RankTwoTensor sigma_s =
      (_lambda_const[_qp] - _damaged_modulus[_qp] / xi) * I1 * I
    + (2.0 * _shear_modulus[_qp] - _damaged_modulus[_qp] * xi) * Ee;

  RankTwoTensor sigma_b =
      (2.0 * _a2[_qp] + _a1[_qp] / xi + 3.0 * _a3[_qp] * xi) * I1 * I
    + (2.0 * _a0[_qp] + _a1[_qp] * xi - _a3[_qp] * std::pow(xi, 3)) * Ee;

  RankTwoTensor sigma_total = (1.0 - _B_breakagevar[_qp]) * sigma_s + _B_breakagevar[_qp] * sigma_b;

  // plastic PK2-like stress & deviatoric measure
  RankTwoTensor Tp = Fe.transpose() * Fe * sigma_total;
  RankTwoTensor Tau = Tp - 0.3333333333 * Tp.trace() * I;

  // tangent modulus using current Ee/I1/I2
  RankFourTensor tangent;
  computeQpTangentModulus(tangent, I1, I2, xi, Ee);

  // --- save properties ---
  _Fp[_qp] = Fp_updated;
  _Fe[_qp] = Fe;
  _Ee[_qp] = Ee;
  _Ep[_qp] = Ep;
  _E[_qp]  = E;

  _I1[_qp] = I1;
  _I2[_qp] = I2;
  _xi[_qp] = xi;

  _S[_qp]  = sigma_total;
  _Tp[_qp] = Tp;
  _Tau[_qp]= Tau;
  _C[_qp]  = tangent;

  // deformation rate & deviatoric strain-rate scalar
  computeDeviatroicStrainRateTensor();
}

RankTwoTensor
ComputeLagrangianDamageBreakageStressPK2DiffusedDebug::computeQpFp()
{
  // Apply power operation on every element of Tau (m2 assumed = 1 here)
  RankTwoTensor Tau_old_power_m2 = _Tau_old[_qp];

  // Plastic deformation rate at t_{n+1} using t_{n} quantities
  RankTwoTensor Dp = _C_g[_qp] * std::pow(_B_breakagevar_old[_qp], _m1[_qp]) * Tau_old_power_m2;

  // Cp = I - Dp dt
  RankTwoTensor Cp = RankTwoTensor::Identity() - Dp * _dt;

  // Implicit Euler update
  RankTwoTensor Fp_updated = Cp.inverse() * _Fp_old[_qp];

  // Save plastic deformation rate
  _Dp[_qp] = Dp;

  return Fp_updated;
}

void
ComputeLagrangianDamageBreakageStressPK2DiffusedDebug::computeQpTangentModulus(RankFourTensor & tangent,
                                                                          Real I1,
                                                                          Real I2,
                                                                          Real xi,
                                                                          RankTwoTensor Ee)
{
  // Use consistent values - same as in stress computation
  Real lambda_out         = _lambda_const[_qp];
  Real shear_modulus_out  = _shear_modulus[_qp];
  Real gamma_damaged_out  = _damaged_modulus[_qp];

  Real a0 = _a0[_qp];
  Real a1 = _a1[_qp];
  Real a2 = _a2[_qp];
  Real a3 = _a3[_qp];

  const Real adjusted_I2 = std::max(I2, 1e-12);
  const Real sqrt_I2     = std::sqrt(adjusted_I2);
  const RankTwoTensor identity = RankTwoTensor::Identity();

  // dξ/dE
  RankTwoTensor dxidE_tensor;
  for (unsigned int k = 0; k < 3; ++k)
    for (unsigned int l = 0; l < 3; ++l)
      dxidE_tensor(k, l) = identity(k, l) / sqrt_I2 - I1 * Ee(k, l) / std::pow(adjusted_I2, 1.5);

  // d(1/ξ)/dE = -1/ξ^2 dξ/dE
  const RankTwoTensor dxim1dE_tensor = dxidE_tensor * (-1.0 / (xi * xi));

  const Real lambda_term = lambda_out - gamma_damaged_out / xi;
  const Real shear_term  = 2.0 * shear_modulus_out - gamma_damaged_out * xi;

  RankFourTensor dSsdE; dSsdE.zero();
  for (unsigned int i = 0; i < 3; ++i)
    for (unsigned int j = 0; j < 3; ++j)
      for (unsigned int k = 0; k < 3; ++k)
        for (unsigned int l = 0; l < 3; ++l)
        {
          // (λ - γ/ξ) δ_kl δ_ij
          dSsdE(i, j, k, l) += lambda_term * identity(i, j) * identity(k, l);
          // (-γ d(1/ξ)/dE_kl) I1 δ_ij
          dSsdE(i, j, k, l) -= gamma_damaged_out * dxim1dE_tensor(k, l) * I1 * identity(i, j);
          // (2μ - γξ) ∂E_ij/∂E_kl
          Real I4_ijkl = 0.5 * (identity(i, k) * identity(j, l) + identity(i, l) * identity(j, k));
          dSsdE(i, j, k, l) += shear_term * I4_ijkl;
          // (-γ dξ/dE_kl) E_ij
          dSsdE(i, j, k, l) -= gamma_damaged_out * dxidE_tensor(k, l) * Ee(i, j);
        }

  const Real coeff2_b = 2.0 * a2 + a1 / xi + 3.0 * a3 * xi;
  const Real coeff4_b = 2.0 * a0 + a1 * xi - a3 * xi * xi * xi;

  RankFourTensor dSbdE; dSbdE.zero();
  for (unsigned int i = 0; i < 3; ++i)
    for (unsigned int j = 0; j < 3; ++j)
      for (unsigned int k = 0; k < 3; ++k)
        for (unsigned int l = 0; l < 3; ++l)
        {
          // (2a2 + a1/ξ + 3a3 ξ) δ_kl δ_ij
          dSbdE(i, j, k, l) += coeff2_b * identity(i, j) * identity(k, l);
          // a1 d(1/ξ)/dE_kl I1 δ_ij
          dSbdE(i, j, k, l) += a1 * dxim1dE_tensor(k, l) * I1 * identity(i, j);
          // 3a3 dξ/dE_kl I1 δ_ij
          dSbdE(i, j, k, l) += 3.0 * a3 * dxidE_tensor(k, l) * I1 * identity(i, j);
          // (2a0 + a1 ξ - a3 ξ^3) ∂E_ij/∂E_kl
          Real I4_ijkl = 0.5 * (identity(i, k) * identity(j, l) + identity(i, l) * identity(j, k));
          dSbdE(i, j, k, l) += coeff4_b * I4_ijkl;
          // a1 dξ/dE_kl E_ij
          dSbdE(i, j, k, l) += a1 * dxidE_tensor(k, l) * Ee(i, j);
          // -a3 d(ξ^3)/dE_kl E_ij  with d(ξ^3)/dE = 3 ξ^2 dξ/dE
          dSbdE(i, j, k, l) -= a3 * 3.0 * xi * xi * dxidE_tensor(k, l) * Ee(i, j);
        }

  tangent = dSsdE * (1.0 - _B_breakagevar[_qp]) + dSbdE * _B_breakagevar[_qp];
}

void
ComputeLagrangianDamageBreakageStressPK2DiffusedDebug::computeDmatrix()
{
  _D[_qp] = 0.5 * ( _velgrad_L[_qp] + _velgrad_L[_qp].transpose() );
}

void
ComputeLagrangianDamageBreakageStressPK2DiffusedDebug::computeDeviatroicStrainRateTensor()
{
  //Compute D matrix
  computeDmatrix();

  // E_dot = F^T D F
  RankTwoTensor E_dot = _F[_qp].transpose() * _D[_qp] * _F[_qp];

  // deviatoric part
  RankTwoTensor E_dev_dot = E_dot - (1.0/3.0) * E_dot.trace() * RankTwoTensor::Identity();

  // J2_dot = 1/2 * E_dev_dot:E_dev_dot
  Real J2_dot = 0.0;
  for (unsigned int i = 0; i < 3; ++i)
    for (unsigned int j = 0; j < 3; ++j)
      J2_dot += 0.5 * E_dev_dot(i,j) * E_dev_dot(i,j);

  _deviatroic_strain_rate[_qp] = std::sqrt(2.0/3.0 * J2_dot);
}
