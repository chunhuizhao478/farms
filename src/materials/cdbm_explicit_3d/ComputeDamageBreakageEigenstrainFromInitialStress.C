//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ComputeDamageBreakageEigenstrainFromInitialStress.h"
#include "RankTwoTensor.h"
#include "Function.h"
#include "Conversion.h" // for stringify

#include <algorithm>
#include <cmath>

registerMooseObject("farmsApp", ComputeDamageBreakageEigenstrainFromInitialStress);

InputParameters
ComputeDamageBreakageEigenstrainFromInitialStress::validParams()
{
  InputParameters params = ComputeEigenstrainBase::validParams();
  params.addClassDescription("Computes an eigenstrain from an initial stress");
  params.addRequiredParam<std::vector<FunctionName>>(
      "initial_stress",
      "A list of functions describing the initial Cauchy stress. There must be 9 of these, "
      "ordered as xx, xy, xz, yx, yy, yz, zx, zy, zz. To compute the eigenstrain correctly, "
      "your elasticity tensor should not be time-varying in the first timestep.");
  params.addCoupledVar("initial_stress_aux",
                       "A list of 9 AuxVariables describing the initial stress. If provided, "
                       "each of these is multiplied by its corresponding initial_stress function "
                       "to obtain the relevant component of initial stress.");
  params.addParam<std::string>("base_name",
                               "The base_name for the elasticity tensor that will be used to "
                               "compute strain from stress. Do not provide any base_name if your "
                               "elasticity tensor does not use one.");
  params.addRequiredParam<Real>("lambda_o", "initial lambda value");
  params.addRequiredParam<Real>("shear_modulus_o", "initial shear modulus value");
  params.addRequiredParam<Real>("xi_o", "xi_o value");
  return params;
}

ComputeDamageBreakageEigenstrainFromInitialStress::ComputeDamageBreakageEigenstrainFromInitialStress(
    const InputParameters & parameters)
  : ComputeEigenstrainBase(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _eigenstrain_old(getMaterialPropertyOld<RankTwoTensor>(_eigenstrain_name)),
    _ini_aux_provided(isParamValid("initial_stress_aux")),
    _ini_aux(_ini_aux_provided ? coupledValues("initial_stress_aux")
                               : std::vector<const VariableValue *>()),
    _lambda_o(getParam<Real>("lambda_o")),
    _shear_modulus_o(getParam<Real>("shear_modulus_o")),
    _xi_o(getParam<Real>("xi_o")),
    _initial_damage_val(getMaterialPropertyByName<Real>("initial_damage"))
{
  const std::vector<FunctionName> & fcn_names(
      getParam<std::vector<FunctionName>>("initial_stress"));
  const std::size_t num = fcn_names.size();

  if (num != LIBMESH_DIM * LIBMESH_DIM)
    paramError("initial_stress",
               "ComputeDamageBreakageEigenstrainFromInitialStress: " +
                   Moose::stringify(LIBMESH_DIM * LIBMESH_DIM) +
                   " initial stress functions must be provided. You supplied " +
                   Moose::stringify(num) + "\n");

  _initial_stress_fcn.resize(num);
  for (unsigned i = 0; i < num; ++i)
    _initial_stress_fcn[i] = &getFunctionByName(fcn_names[i]);

  if (_ini_aux_provided)
  {
    const std::size_t aux_size = coupledComponents("initial_stress_aux");
    if (aux_size != LIBMESH_DIM * LIBMESH_DIM)
      paramError("initial_stress_aux",
                 "ComputeDamageBreakageEigenstrainFromInitialStress: If you supply "
                 "initial_stress_aux, " +
                     Moose::stringify(LIBMESH_DIM * LIBMESH_DIM) +
                     " values must be given. You supplied " + Moose::stringify(aux_size) + "\n");
  }
}

void
ComputeDamageBreakageEigenstrainFromInitialStress::computeQpEigenstrain()
{
  // Only initialize at first step – afterwards carry old value
  if (_t_step != 1)
  {
    _eigenstrain[_qp] = _eigenstrain_old[_qp];
    return;
  }

  //-----------------------------------------------------------
  // 1) Read initial stress σ⁰ and symmetrize
  //-----------------------------------------------------------
  RankTwoTensor sigma;
  for (unsigned i = 0; i < LIBMESH_DIM; ++i)
    for (unsigned j = 0; j < LIBMESH_DIM; ++j)
    {
      // index order: i*dim + j -> xx, xy, xz, yx, yy, yz, zx, zy, zz
      sigma(i, j) = _initial_stress_fcn[i * LIBMESH_DIM + j]->value(_t, _q_point[_qp]);
      if (_ini_aux_provided)
        sigma(i, j) *= (*_ini_aux[i * LIBMESH_DIM + j])[_qp];
    }

  RankTwoTensor sigma_sym = sigma;
  sigma_sym += sigma.transpose();
  sigma_sym *= 0.5;

  //-----------------------------------------------------------
  // 2) Damage-dependent parameters (B=0 case)
  //-----------------------------------------------------------
  Real gamma_r = computegammar();
  if (!std::isfinite(gamma_r))
  {
    mooseWarning("computegammar() returned non-finite at qp=", _qp, " -> using 0");
    gamma_r = 0.0;
  }

  const Real damage = _initial_damage_val[_qp]; // d
  const Real gamma  = damage * gamma_r;         // γ = d γ_r
  const Real mu     = _shear_modulus_o + _xi_o * damage * gamma_r; // μ(d)
  const Real lambda = _lambda_o;                // keep λ constant here

  //-----------------------------------------------------------
  // 3) Newton iteration for ξ (with sign-preserving floor & bounds)
  //-----------------------------------------------------------
  const Real S1 = sigma_sym.trace();                        // σ_kk
  const Real S2 = sigma_sym.doubleContraction(sigma_sym);   // σ:σ
  const Real n  = static_cast<Real>(LIBMESH_DIM);

  auto keep_away_from_zero = [](Real x, Real eps)
  {
    if (std::abs(x) >= eps)
      return x;
    return (x < 0.0) ? -eps : eps; // preserve sign; if x==0 choose +eps
  };

  const Real xi_eps = 1e-12;
  const Real xi_max = std::sqrt(3.0) - 1e-10; // physical |ξ| ≤ √3

  // initial guess: sign of S1 retained if possible
  Real xi = (std::abs(S1) > 1e-16 && std::abs(S2) > 1e-16) ? (S1 / std::sqrt(S2)) : _xi_o;
  xi = keep_away_from_zero(xi, xi_eps);
  xi = std::clamp(xi, -xi_max, xi_max);

  // Helper to compute residual F(ξ) = ξ − I1(ξ)/sqrt(I2(ξ))
  auto compute_residual = [&](Real x, Real & F_out) -> bool
  {
    const Real a = lambda - gamma / x;
    const Real b = 2.0 * mu - gamma * x;
    const Real D = n * a + b;

    if (std::abs(b) < 1e-18 || std::abs(D) < 1e-18)
      return false;

    const Real A = 1.0 / b;
    const Real Bcomp = -a / (b * D);

    const Real I1 = S1 / D;
    const Real I2 = A * A * S2 + (2.0 * A * Bcomp + n * Bcomp * Bcomp) * (S1 * S1);
    if (!(I2 > 0.0) || !std::isfinite(I2))
      return false;

    const Real sqrtI2 = std::sqrt(I2);
    F_out = x - I1 / sqrtI2;
    return std::isfinite(F_out);
  };

  // Newton loop
  for (unsigned it = 0; it < 20; ++it)
  {
    // Build at current xi
    const Real a = lambda - gamma / xi;
    const Real b = 2.0 * mu - gamma * xi;
    const Real D = n * a + b;

    if (std::abs(b) < 1e-18 || std::abs(D) < 1e-18)
    {
      mooseWarning("Ill-conditioned operator at qp=", _qp, " (b=", b, ", D=", D, ")");
      break;
    }

    const Real A = 1.0 / b;
    const Real Bcomp = -a / (b * D);

    const Real I1 = S1 / D;
    Real I2 = A * A * S2 + (2.0 * A * Bcomp + n * Bcomp * Bcomp) * (S1 * S1);
    I2 = std::max(I2, 1e-32);
    const Real sqrtI2 = std::sqrt(I2);

    Real F = xi - I1 / sqrtI2;
    if (!std::isfinite(F))
    {
      mooseWarning("Residual non-finite at qp=", _qp, " (xi=", xi, ")");
      break;
    }
    if (std::abs(F) < 1e-12)
      break; // converged

    // Derivative dF/dξ
    const Real da =  gamma / (xi * xi);
    const Real db = -gamma;
    const Real dD = n * da + db;

    const Real dI1 = -S1 * dD / (D * D);

    const Real dA = -db / (b * b); // = gamma / b^2
    const Real dB = (-da * b * D + a * (db * D + b * dD)) / (b * b * D * D);
    const Real dK = 2.0 * (dA * Bcomp + A * dB) + 2.0 * n * Bcomp * dB; // K' where K=2AB+nB^2
    const Real dI2 = 2.0 * A * dA * S2 + dK * (S1 * S1);

    const Real denomI2 = std::max(I2, 1e-32);
    const Real dF = 1.0 - dI1 / sqrtI2 + 0.5 * I1 * dI2 / (denomI2 * sqrtI2);

    if (!std::isfinite(dF) || std::abs(dF) < 1e-18)
    {
      mooseWarning("dF breakdown at qp=", _qp, " (xi=", xi, ")");
      break;
    }

    // Newton step with lightweight backtracking line search
    Real step = -F / dF;
    Real best_xi = xi;
    Real best_F  = F;

    // Try the full step first
    Real xi_try = std::clamp(keep_away_from_zero(xi + step, xi_eps), -xi_max, xi_max);
    Real F_try;
    if (compute_residual(xi_try, F_try) && std::abs(F_try) < std::abs(best_F))
    {
      best_xi = xi_try;
      best_F  = F_try;
    }
    else
    {
      // Backtrack up to 8 times
      Real t = 0.5;
      for (unsigned ls = 0; ls < 8; ++ls)
      {
        xi_try = std::clamp(keep_away_from_zero(xi + t * step, xi_eps), -xi_max, xi_max);
        if (compute_residual(xi_try, F_try) && std::abs(F_try) < std::abs(best_F))
        { best_xi = xi_try; best_F = F_try; break; }
        t *= 0.5;
      }
    }

    xi = best_xi;
    if (std::abs(best_F) < 1e-12)
      break;
  }

  //-----------------------------------------------------------
  // 4) Compute compliance constants and eigenstrain ε^e
  //-----------------------------------------------------------
  const Real a = lambda - gamma / xi;
  const Real b = 2.0 * mu - gamma * xi;
  const Real D = n * a + b;

  if (std::abs(b) < 1e-18 || std::abs(D) < 1e-18)
  {
    mooseWarning("Ill-conditioned inversion at qp=", _qp, " (b=", b, ", D=", D, "). "
                 "Setting eigenstrain to zero.");
    _eigenstrain[_qp].zero();
    return;
  }

  const Real A = 1.0 / b;
  const Real Bcomp = -a / (b * D);

  if (!std::isfinite(A) || !std::isfinite(Bcomp))
  {
    mooseWarning("Non-finite compliance at qp=", _qp, " (A=", A, ", Bcomp=", Bcomp, "). "
                 "Setting eigenstrain to zero.");
    _eigenstrain[_qp].zero();
    return;
  }

  RankTwoTensor eps = sigma_sym;                 // start with sym σ_ij
  eps *= A;                                      // A σ_ij
  eps += Bcomp * S1 * RankTwoTensor::Identity(); // + Bcomp tr(σ) δ_ij

  if (!(std::isfinite(eps(0,0)) && std::isfinite(eps(1,1)) && std::isfinite(eps(2,2)) &&
        std::isfinite(eps(0,1)) && std::isfinite(eps(0,2)) && std::isfinite(eps(1,2))))
  {
    mooseWarning("Computed strain has non-finite entries at qp=", _qp, ". Setting to zero.");
    _eigenstrain[_qp].zero();
    return;
  }

  // Store as negative eigenstrain (initial strain we want to cancel)
  _eigenstrain[_qp] = -eps;
}

Real
ComputeDamageBreakageEigenstrainFromInitialStress::computegammar()
{
  // Calculate each part with guards and branch selection
  Real denom = 2.0 * (std::pow(_xi_o, 2) - 3.0);
  if (std::abs(denom) < 1e-18)
  {
    mooseWarning("computegammar(): denominator near zero for xi_o=", _xi_o, ". Returning 0.");
    return 0.0;
  }

  Real term1 = -_xi_o * (-_lambda_o * std::pow(_xi_o, 2) + 6.0 * _lambda_o + 2.0 * _shear_modulus_o);

  Real r1 = (_lambda_o * std::pow(_xi_o, 2) + 2.0 * _shear_modulus_o);
  Real r2 = (_lambda_o * std::pow(_xi_o, 4) - 12.0 * _lambda_o * std::pow(_xi_o, 2) + 36.0 * _lambda_o
             - 6.0 * _shear_modulus_o * std::pow(_xi_o, 2) + 24.0 * _shear_modulus_o);
  Real rad = r1 * r2;

  if (rad < 0.0)
  {
    mooseWarning("computegammar(): negative radicand (", rad, "). Clamping to zero.");
    rad = 0.0;
  }

  Real term2_sqrt = std::sqrt(rad);

  Real gamma_r_minus = (term1 - term2_sqrt) / denom;
  Real gamma_r_plus  = (term1 + term2_sqrt) / denom;

  // Prefer the branch that keeps b0 = 2μ0 − γ_r ξ0 positive
  Real b_minus = 2.0 * _shear_modulus_o - gamma_r_minus * _xi_o;
  if (b_minus > 0.0)
    return gamma_r_minus;
  else
    return gamma_r_plus;
}
