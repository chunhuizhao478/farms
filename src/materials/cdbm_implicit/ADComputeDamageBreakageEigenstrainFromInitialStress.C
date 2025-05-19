//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADComputeDamageBreakageEigenstrainFromInitialStress.h"
#include "RankTwoTensor.h"
#include "Function.h"
#include "Conversion.h" // for stringify

registerMooseObject("farmsApp", ADComputeDamageBreakageEigenstrainFromInitialStress);

InputParameters
ADComputeDamageBreakageEigenstrainFromInitialStress::validParams()
{
  InputParameters params = ADComputeEigenstrainBase::validParams();
  params.addClassDescription("Computes an eigenstrain from an initial strain");
  params.addRequiredParam<std::vector<FunctionName>>(
      "initial_stress",
      "A list of functions describing the initial stress.  There must be 9 of these, corresponding "
      "to the xx, yx, zx, xy, yy, zy, xz, yz, zz components respectively.  To compute the "
      "eigenstrain correctly, your elasticity tensor should not be time-varying in the first "
      "timestep");
  params.addCoupledVar("initial_stress_aux",
                       "A list of 9 AuxVariables describing the initial stress.  If provided, each "
                       "of these is multiplied by its corresponding initial_stress function to "
                       "obtain the relevant component of initial strain.");
  params.addParam<std::string>("base_name",
                               "The base_name for the elasticity tensor that will be "
                               "used to compute strain from stress.  Do not provide "
                               "any base_name if your elasticity tensor does not use "
                               "one.");
  params.addRequiredParam<Real>("lambda_o","initial lambda value");
  params.addRequiredParam<Real>("shear_modulus_o","initial shear modulus value");
  params.addRequiredParam<Real>("xi_o","xi_o value");
  return params;
}

ADComputeDamageBreakageEigenstrainFromInitialStress::ADComputeDamageBreakageEigenstrainFromInitialStress(
    const InputParameters & parameters)
  : ADComputeEigenstrainBase(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _eigenstrain_old(getMaterialPropertyOld<RankTwoTensor>(_eigenstrain_name)),
    _ini_aux_provided(isParamValid("initial_stress_aux")),
    _ini_aux(_ini_aux_provided ? coupledValues("initial_stress_aux")
                               : std::vector<const VariableValue *>()),
    _lambda_o(getParam<Real>("lambda_o")),
    _shear_modulus_o(getParam<Real>("shear_modulus_o")),
    _xi_o(getParam<Real>("xi_o")),
    _initial_damage_val(getADMaterialPropertyByName<Real>("initial_damage"))
{
  const std::vector<FunctionName> & fcn_names(
      getParam<std::vector<FunctionName>>("initial_stress"));
  const std::size_t num = fcn_names.size();

  if (num != LIBMESH_DIM * LIBMESH_DIM)
    paramError(
        "initial_stress",
        "ADComputeDamageBreakageEigenstrainFromInitialStress: " + Moose::stringify(LIBMESH_DIM * LIBMESH_DIM) +
            " initial strain functions must be provided.  You supplied " + Moose::stringify(num) +
            "\n");

  _initial_stress_fcn.resize(num);
  for (unsigned i = 0; i < num; ++i)
    _initial_stress_fcn[i] = &getFunctionByName(fcn_names[i]);

  if (_ini_aux_provided)
  {
    const std::size_t aux_size = coupledComponents("initial_stress_aux");
    if (aux_size != LIBMESH_DIM * LIBMESH_DIM)
      paramError("initial_stress_aux",
                 "ADComputeDamageBreakageEigenstrainFromInitialStress: If you supply initial_stress_aux, " +
                     Moose::stringify(LIBMESH_DIM * LIBMESH_DIM) +
                     " values must be given.  You supplied " + Moose::stringify(aux_size) + "\n");
  }
}

void
ADComputeDamageBreakageEigenstrainFromInitialStress::computeQpEigenstrain()
{
  // only initialise at the first time step – afterwards carry old value
  if (_t_step != 1)
  {
    _eigenstrain[_qp] = _eigenstrain_old[_qp];
    return;
  }

  //-----------------------------------------------------------
  // 1. Read initial stress σ⁰
  //-----------------------------------------------------------
  ADRankTwoTensor sigma; // σ_ij
  for (unsigned i = 0; i < LIBMESH_DIM; ++i)
    for (unsigned j = 0; j < LIBMESH_DIM; ++j)
    {
      sigma(i, j) = _initial_stress_fcn[i * LIBMESH_DIM + j]->value(_t, _q_point[_qp]);
      if (_ini_aux_provided)
        sigma(i, j) *= (*_ini_aux[i * LIBMESH_DIM + j])[_qp];
    }

  //-----------------------------------------------------------
  // 2. Damage‑dependent parameters
  //-----------------------------------------------------------
  ADReal gamma_r = computegammar();                                    // γ_r
  ADReal damage  = _initial_damage_val[_qp];                           // d
  ADReal gamma   = damage * gamma_r;                                   // γ = d γ_r
  ADReal mu      = _shear_modulus_o + _xi_o * damage * gamma_r;        // μ(ξ₀,d)
  ADReal lambda  = _lambda_o;                                          // λ (kept constant here)

  //-----------------------------------------------------------
  // 3. Newton iteration for ξ
  //-----------------------------------------------------------
  const ADReal S1 = sigma.trace();                 // σ_kk
  const ADReal S2 = (sigma * sigma).trace();       // σ_ij σ_ij

  // initial guess: use reference value or ratio of invariants if possible
  ADReal xi = (std::abs(S1) > 1e-16 && std::abs(S2) > 1e-16) ? S1 / std::sqrt(S2) : _xi_o;

  for (unsigned it = 0; it < 10; ++it)
  {
    // a = λ - γ/ξ,  b = 2μ - γ ξ
    ADReal a = lambda - gamma / xi;
    ADReal b = 2.0 * mu - gamma * xi;
    ADReal denom = 3.0 * a + b; // D

    // Invariants I1, I2 as functions of ξ
    ADReal I1 = S1 / denom;

    // coefficient C(ξ) used in I2 expression
    ADReal C = (2.0 * a) / (b * b * denom) + (3.0 * a * a) / (b * b * denom * denom);
    ADReal I2 = S2 / (b * b) + C * S1 * S1;

    ADReal sqrtI2 = std::sqrt(I2 + 1e-32);

    // Residual F(ξ)
    ADReal F = xi - I1 / sqrtI2;

    if (std::abs(MetaPhysicL::raw_value(F)) < 1e-12) // converged (raw_value for ADReal tolerance)
      break;

    //-------------------------------------------------------
    // Derivative F'(ξ)
    //-------------------------------------------------------
    ADReal da =  gamma / (xi * xi);   // a'
    ADReal db = -gamma;               // b'
    ADReal dDenom = 3.0 * da + db;    // D'
    ADReal dI1 = -S1 * dDenom / (denom * denom);

    // derivative of C(ξ)
    ADReal dC = (2.0 * da) / (b * b * denom)
               - (4.0 * a * db) / (b * b * b * denom)
               - (2.0 * a * dDenom) / (b * b * denom * denom)
               + (6.0 * a * da) / (b * b * denom * denom)
               - (6.0 * a * a * db) / (b * b * b * denom * denom)
               - (6.0 * a * a * dDenom) / (b * b * denom * denom * denom);

    ADReal dI2 = -2.0 * S2 * db / (b * b * b) + dC * S1 * S1;

    ADReal dF = 1.0 - dI1 / sqrtI2 + 0.5 * I1 * dI2 / (I2 * sqrtI2);

    xi -= F / dF; // Newton update
  }

  if (it == 10)
    mooseError("ADComputeDamageBreakageEigenstrainFromInitialStress: "
               "Newton iteration for xi did not converge.  "
               "xi = " + Moose::stringify(xi) + "\n");

  //-----------------------------------------------------------
  // 4. Compute compliance constants A,B and eigenstrain ε^e
  //-----------------------------------------------------------
  ADReal a = lambda - gamma / xi;
  ADReal b = 2.0 * mu - gamma * xi;
  ADReal A = 1.0 / b;
  ADReal B = -a / (b * (3.0 * a + b));

  ADRankTwoTensor eps = sigma;     // start with σ_ij
  eps *= A;                        // A σ_ij
  eps += B * S1 * RankTwoTensor::Identity(); // + B tr(σ) δ_ij

  // Store as negative eigen‑strain (initial strain we want to cancel)
  _eigenstrain[_qp] = -eps;
}

ADReal
ADComputeDamageBreakageEigenstrainFromInitialStress::computegammar()
{
  // Calculate each part of the expression
  ADReal term1 = -_xi_o * (-_lambda_o * std::pow(_xi_o, 2) + 6 * _lambda_o + 2 * _shear_modulus_o);
  ADReal term2_sqrt = std::sqrt((_lambda_o * std::pow(_xi_o, 2) + 2 * _shear_modulus_o) * 
                            (_lambda_o * std::pow(_xi_o, 4) - 12 * _lambda_o * std::pow(_xi_o, 2) + 36 * _lambda_o 
                            - 6 * _shear_modulus_o * std::pow(_xi_o, 2) + 24 * _shear_modulus_o));
  ADReal denominator = 2 * (std::pow(_xi_o, 2) - 3);
  
  // Calculate gamma_r
  ADReal gamma_r = (term1 - term2_sqrt) / denominator;
  
  return gamma_r;
}
