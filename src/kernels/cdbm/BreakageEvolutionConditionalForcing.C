//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "BreakageEvolutionConditionalForcing.h"

registerMooseObject("farmsApp", BreakageEvolutionConditionalForcing);

InputParameters
BreakageEvolutionConditionalForcing::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("This kernel implements the damage evolution forcing term: $(1-B) psi f(alpha, xi)$");
  params.addRequiredCoupledVar("coupled", "The damage variable");
  return params;
}

BreakageEvolutionConditionalForcing::BreakageEvolutionConditionalForcing(const InputParameters & parameters)
  : Kernel(parameters),
    _alpha_var(coupled("coupled")),
    _alpha(coupledValue("coupled")),
    _Cd(getMaterialProperty<Real>("Cd")),
    _CdCb_multiplier(getMaterialProperty<Real>("CdCb_multiplier")),
    _I2(getMaterialProperty<Real>("I2")),
    _xi(getMaterialProperty<Real>("xi")),
    _xi_0(getMaterialProperty<Real>("xi_0")),
    _xi_1(getMaterialProperty<Real>("xi_1")),
    _xi_d(getMaterialProperty<Real>("xi_d")),
    _xi_min(getMaterialProperty<Real>("xi_min")),
    _xi_max(getMaterialProperty<Real>("xi_max")),
    _C1(getMaterialProperty<Real>("C_1")),
    _C2(getMaterialProperty<Real>("C_2")),
    _lambda_o(getMaterialProperty<Real>("lambda_o")),
    _shear_modulus_o(getMaterialProperty<Real>("shear_modulus_o")),
    _CBH_constant(getMaterialProperty<Real>("CBH_constant")),
    _beta_width(getMaterialProperty<Real>("beta_width")),
    _gamma_damaged_r(getMaterialProperty<Real>("gamma_damaged_r"))
{
}

Real
BreakageEvolutionConditionalForcing::computeQpResidual()
{
  return -1.0 * _test[_i][_qp] * computebreakageevolutionforcingterm();
}

Real
BreakageEvolutionConditionalForcing::computeQpJacobian()
{
  // The partial derivative of the integrand wrt _u:
  const Real dforcing_dB = computebreakageevolutionforcingterm_derivative();

  // In MOOSE, the standard Jacobian for a Kernel is:
  //   J = (dF/dalpha) * test[_i][_qp] * phi[_j][_qp]
  // Here we also have the factor -1.0 from your residual:
  return -1.0 * _test[_i][_qp] * dforcing_dB * _phi[_j][_qp];
}

Real
BreakageEvolutionConditionalForcing::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _alpha_var)
    return -1.0 * _test[_i][_qp] * computebreakageevolutionforcingterm_derivative_walpha();
  else
    return 0.0;
}

Real
BreakageEvolutionConditionalForcing::computebreakageevolutionforcingterm()
{
  /* compute C_B based on C_d */
  Real C_B = _CdCb_multiplier[_qp] * _Cd[_qp];

  //alphacr function
  Real alphacr;
  if ( _xi[_qp] < _xi_0[_qp] && _xi[_qp] >= _xi_min[_qp] ){ alphacr = 1.0;} 
  else if ( _xi[_qp] > _xi_0[_qp] && _xi[_qp] <= _xi_1[_qp] ){ alphacr = alphacr_root1(_xi[_qp]);}
  else if ( _xi[_qp] > _xi_1[_qp] && _xi[_qp] <= _xi_max[_qp] ){ alphacr = alphacr_root2(_xi[_qp]);}
  else{std::cout<<"xi: "<<_xi[_qp]<<std::endl; mooseError("xi exceeds the maximum allowable range!");}

  //compute forcing func
  Real Prob = 1.0 / ( std::exp( (alphacr - _alpha[_qp]) / _beta_width[_qp] ) + 1.0 );

  Real B_forcingterm;
  if (_xi[_qp] >= _xi_d[_qp]){
    B_forcingterm = 1.0 * C_B * Prob * (1-_u[_qp]) * _I2[_qp] * (_xi[_qp] - _xi_d[_qp]); 
  }
  else if (_xi[_qp] < _xi_d[_qp]){
    B_forcingterm = 1.0 * _CBH_constant[_qp] * _I2[_qp] * (_xi[_qp] - _xi_d[_qp]);
  }
  else{
    mooseError("xi is OUT-OF-RANGE!.");   
  }

  return B_forcingterm;
}

Real
BreakageEvolutionConditionalForcing::computebreakageevolutionforcingterm_derivative()
{
  /* compute C_B based on C_d */
  Real C_B = _CdCb_multiplier[_qp] * _Cd[_qp];

  //alphacr function
  Real alphacr;
  if ( _xi[_qp] < _xi_0[_qp] && _xi[_qp] >= _xi_min[_qp] ){ alphacr = 1.0;} 
  else if ( _xi[_qp] > _xi_0[_qp] && _xi[_qp] <= _xi_1[_qp] ){ alphacr = alphacr_root1(_xi[_qp]);}
  else if ( _xi[_qp] > _xi_1[_qp] && _xi[_qp] <= _xi_max[_qp] ){ alphacr = alphacr_root2(_xi[_qp]);}
  else{std::cout<<"xi: "<<_xi[_qp]<<std::endl; mooseError("xi exceeds the maximum allowable range!");}

  //compute forcing func
  Real Prob = 1.0 / ( std::exp( (alphacr - _alpha[_qp]) / _beta_width[_qp] ) + 1.0 );

  Real dforcing_dB;
  if (_xi[_qp] >= _xi_d[_qp]){
    dforcing_dB = -1.0 * C_B * Prob * _I2[_qp] * (_xi[_qp] - _xi_d[_qp]); 
  }
  else if (_xi[_qp] < _xi_d[_qp]){
    dforcing_dB = 0.0;
  }
  else{
    mooseError("xi is OUT-OF-RANGE!.");   
  }

  return dforcing_dB;
}

Real
BreakageEvolutionConditionalForcing::computebreakageevolutionforcingterm_derivative_walpha()
{
  /* compute C_B based on C_d */
  Real C_B = _CdCb_multiplier[_qp] * _Cd[_qp];

  //alphacr function
  Real alphacr;
  if ( _xi[_qp] < _xi_0[_qp] && _xi[_qp] >= _xi_min[_qp] ){ alphacr = 1.0;} 
  else if ( _xi[_qp] > _xi_0[_qp] && _xi[_qp] <= _xi_1[_qp] ){ alphacr = alphacr_root1(_xi[_qp]);}
  else if ( _xi[_qp] > _xi_1[_qp] && _xi[_qp] <= _xi_max[_qp] ){ alphacr = alphacr_root2(_xi[_qp]);}
  else{std::cout<<"xi: "<<_xi[_qp]<<std::endl; mooseError("xi exceeds the maximum allowable range!");}

  // Here we assume that the derivative of the logistic function is the only α-dependence.
  // The derivative of 1/(exp((alphacr-α)/β)+1) with respect to α is
  //   (exp((alphacr-α)/β))/(β*(exp((alphacr-α)/β)+1)^2)
  // Multiplying by the other constant factors (and noting the chain rule),
  // a simplified formulation is given here (adjust as needed for your model):
  Real dProb_dalpha = std::exp((alphacr - _alpha[_qp]) / _beta_width[_qp]) /
  (_beta_width[_qp] * std::pow(std::exp((alphacr - _alpha[_qp]) / _beta_width[_qp]) + 1.0, 2));

  Real dforcing_dalpha;
  if (_xi[_qp] >= _xi_d[_qp]){
    dforcing_dalpha = C_B * (1-_u[_qp]) * dProb_dalpha * _I2[_qp] * (_xi[_qp] - _xi_d[_qp]); 
  }
  else if (_xi[_qp] < _xi_d[_qp]){
    dforcing_dalpha = 0.0;
  }
  else{
    mooseError("xi is OUT-OF-RANGE!.");   
  }

  return dforcing_dalpha;
}

// Function for alpha_func_root1
Real 
BreakageEvolutionConditionalForcing::alphacr_root1(Real xi) {

  Real term1 = _lambda_o[_qp] * pow(xi, 3) - 6 * _lambda_o[_qp] * _xi_0[_qp] + 6 * _shear_modulus_o[_qp] * xi - 8 * _shear_modulus_o[_qp] * _xi_0[_qp];
  Real term2 = std::sqrt(_lambda_o[_qp] * _lambda_o[_qp] * pow(xi, 6) 
                            - 12 * _lambda_o[_qp] * _lambda_o[_qp] * pow(xi, 3) * _xi_0[_qp] 
                            + 36 * _lambda_o[_qp] * _lambda_o[_qp] * _xi_0[_qp] * _xi_0[_qp] 
                            + 12 * _lambda_o[_qp] * _shear_modulus_o[_qp] * pow(xi, 4) 
                            - 16 * _lambda_o[_qp] * _shear_modulus_o[_qp] * pow(xi, 3) * _xi_0[_qp] 
                            - 72 * _lambda_o[_qp] * _shear_modulus_o[_qp] * pow(xi, 2) 
                            + 72 * _lambda_o[_qp] * _shear_modulus_o[_qp] * xi * _xi_0[_qp] 
                            + 72 * _lambda_o[_qp] * _shear_modulus_o[_qp] 
                            - 12 * _shear_modulus_o[_qp] * _shear_modulus_o[_qp] * pow(xi, 2) 
                            + 48 * _shear_modulus_o[_qp] * _shear_modulus_o[_qp]);
  Real denominator = 2 * _gamma_damaged_r[_qp] * (3 * pow(xi, 2) - 6 * xi * _xi_0[_qp] + 4 * _xi_0[_qp] * _xi_0[_qp] - 3);
  return (term1 - term2) / denominator;
}

// Function for alpha_func_root2
Real 
BreakageEvolutionConditionalForcing::alphacr_root2(Real xi) {
  return 2 * _shear_modulus_o[_qp] / (_gamma_damaged_r[_qp] * (xi - 2 * _xi_0[_qp]));
}