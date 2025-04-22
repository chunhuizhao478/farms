//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DamageEvolutionConditionalForcing.h"

registerMooseObject("farmsApp", DamageEvolutionConditionalForcing);

InputParameters
DamageEvolutionConditionalForcing::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("This kernel implements the breakage evolution forcing term: $(1-B) psi f(alpha, xi)$");
  params.addRequiredCoupledVar("coupled", "The breakage variable");
  return params;
}

DamageEvolutionConditionalForcing::DamageEvolutionConditionalForcing(const InputParameters & parameters)
  : Kernel(parameters),
    _B_var(coupled("coupled")),
    _B(coupledValue("coupled")),
    _Cd(getMaterialProperty<Real>("Cd")),
    _I2(getMaterialProperty<Real>("I2")),
    _xi(getMaterialProperty<Real>("xi")),
    _xi_0(getMaterialProperty<Real>("xi_0")),
    _C1(getMaterialProperty<Real>("C_1")),
    _C2(getMaterialProperty<Real>("C_2"))
{
}

Real
DamageEvolutionConditionalForcing::computeQpResidual()
{
  return -1.0 * _test[_i][_qp] * (1 - _B[_qp]) * computedamageevolutionforcingterm();
}

Real
DamageEvolutionConditionalForcing::computeQpJacobian()
{
  // The partial derivative of the integrand wrt _u:
  const Real dforcing_dalpha = computedamageevolutionforcingterm_derivative();

  // In MOOSE, the standard Jacobian for a Kernel is:
  //   J = (dF/dalpha) * test[_i][_qp] * phi[_j][_qp]
  // Here we also have the factor -1.0 from your residual:
  return -_test[_i][_qp] * dforcing_dalpha * _phi[_j][_qp];
}

Real
DamageEvolutionConditionalForcing::computeQpOffDiagJacobian(unsigned int jvar)
{
  // This function computes the derivative of the residual with respect to the coupled variable B.
  // Recall that: R = -test * (1-B) * f(u, xi)
  // So, dR/dB = -test * (-f(u, xi)) = test * f(u, xi)
  if (jvar == _B_var)
    return _test[_i][_qp] * computedamageevolutionforcingterm();
  else
    return 0.0;
}

Real
DamageEvolutionConditionalForcing::computedamageevolutionforcingterm()
{
  Real alpha_forcingterm;
  if (_xi[_qp] >= _xi_0[_qp]){
    alpha_forcingterm = ( _Cd[_qp] * _I2[_qp] * ( _xi[_qp] - _xi_0[_qp] ) );
  }
  else if (_xi[_qp] < _xi_0[_qp]){
    //_u[_qp] == alpha
    alpha_forcingterm = ( _C1[_qp] * std::exp(_u[_qp]/_C2[_qp]) * _I2[_qp] * ( _xi[_qp] - _xi_0[_qp] ) );
  }
  else{
    mooseError("xi is OUT-OF-RANGE!.");   
  }

  return alpha_forcingterm;
}

Real
DamageEvolutionConditionalForcing::computedamageevolutionforcingterm_derivative()
{
  // Derivative of the quantity used in 'computedamageevolutionforcingterm()'
  // with respect to the primary variable alpha = _u[_qp]

  Real dforcing_dalpha = 0.0;

  if (_xi[_qp] >= _xi_0[_qp])
  {
    // For xi >= xi_0, there's no direct dependence on alpha (u), so derivative = 0
    dforcing_dalpha = 0.0;
  }
  else if (_xi[_qp] < _xi_0[_qp])
  {
    // For xi < xi_0, the forcing term is:
    //    (1 - B) * (C1 * exp(u / C2) * I2 * (xi - xi_0))
    // so its derivative w.r.t. u is:
    //
    //    (1 - B) * [C1 / C2 * exp(u / C2) * I2 * (xi - xi_0)]
    //
    dforcing_dalpha =
      (1.0 - _B[_qp]) *
      ( _C1[_qp] * std::exp(_u[_qp] / _C2[_qp]) * _I2[_qp] / _C2[_qp] * (_xi[_qp] - _xi_0[_qp]) );
  }
  else
  {
    // xi == xi_0 is typically a boundary case, pick whichever side you need or set 0
    dforcing_dalpha = 0.0;
  }

  return dforcing_dalpha;
}