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
  if (_t_step == 1) //load the solution from initial
  {
    ADRankTwoTensor initial_stress;
    for (unsigned i = 0; i < LIBMESH_DIM; ++i)
      for (unsigned j = 0; j < LIBMESH_DIM; ++j)
      {
        initial_stress(i, j) = _initial_stress_fcn[i * LIBMESH_DIM + j]->value(_t, _q_point[_qp]);
        if (_ini_aux_provided)
          initial_stress(i, j) *= (*_ini_aux[i * LIBMESH_DIM + j])[_qp];
      }

    // Compute gammar
    ADReal gamma_r = computegammar();

    // Evaluate shear modulus
    ADReal shear_modulus = _shear_modulus_o + _xi_o * _initial_damage_val[_qp] * gamma_r;
    ADReal gamma_damaged_out = _initial_damage_val[_qp] * gamma_r;

    // Compute strain components
    ADReal epsxx = 1.0/(2.0*shear_modulus) * initial_stress(0,0) - _lambda_o / ((2.0 * shear_modulus)*(3*_lambda_o + 2*shear_modulus)) * (initial_stress(0,0) + initial_stress(1,1) + initial_stress(2,2));
    ADReal epsyy = 1.0/(2.0*shear_modulus) * initial_stress(1,1) - _lambda_o / ((2.0 * shear_modulus)*(3*_lambda_o + 2*shear_modulus)) * (initial_stress(0,0) + initial_stress(1,1) + initial_stress(2,2));
    ADReal epszz = 1.0/(2.0*shear_modulus) * initial_stress(2,2) - _lambda_o / ((2.0 * shear_modulus)*(3*_lambda_o + 2*shear_modulus)) * (initial_stress(0,0) + initial_stress(1,1) + initial_stress(2,2));
    ADReal epsxy = 1.0/(2.0*shear_modulus) * initial_stress(0,1);
    ADReal epsxz = 1.0/(2.0*shear_modulus) * initial_stress(0,2);
    ADReal epsyz = 1.0/(2.0*shear_modulus) * initial_stress(1,2);
    
    //invSymm only works for non-AD
    _eigenstrain[_qp](0,0) = epsxx;
    _eigenstrain[_qp](1,1) = epsyy;
    _eigenstrain[_qp](2,2) = epszz;
    _eigenstrain[_qp](0,1) = epsxy; _eigenstrain[_qp](1,0) = epsxy;
    _eigenstrain[_qp](0,2) = epsxz; _eigenstrain[_qp](2,0) = epsxz;
    _eigenstrain[_qp](1,2) = epsyz; _eigenstrain[_qp](2,1) = epsyz;
  }
  else
    _eigenstrain[_qp] = _eigenstrain_old[_qp];
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
