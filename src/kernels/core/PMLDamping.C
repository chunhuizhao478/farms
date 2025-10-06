//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

/*
 * Perfectly Matched Layer (PML) damping kernel for absorbing boundaries
 * Implements damping term for time-domain PML in elastic wave propagation
 */

#include "PMLDamping.h"

registerMooseObject("farmsApp", PMLDamping);

InputParameters
PMLDamping::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("PML damping kernel that adds spatially-varying damping "
                             "in the absorbing layer to eliminate boundary reflections. "
                             "Residual: -d(x,y) * rho * du/dt");
  return params;
}

PMLDamping::PMLDamping(const InputParameters & parameters)
  : Kernel(parameters),
    _u_dot(dot()),
    _density(getMaterialProperty<Real>("density")),
    _pml_damping_coeff(getMaterialProperty<Real>("pml_damping_coeff"))
{
}

Real
PMLDamping::computeQpResidual()
{
  // PML damping term: -d * rho * du/dt
  // This absorbs waves in the PML layer
  return _test[_i][_qp] * _pml_damping_coeff[_qp] * _density[_qp] * _u_dot[_qp];
}

Real
PMLDamping::computeQpJacobian()
{
  // Jacobian w.r.t. u: d(Residual)/du = d * rho * du_dot/du
  // For time integration: du_dot/du = a (depends on time integrator)
  // For explicit central difference: contribution is zero (mass matrix only)
  // For implicit: would need du_dot/du from time integrator

  // Since this is for explicit dynamics with central difference,
  // the Jacobian contribution is zero (diagonal mass matrix is constant)
  return 0.0;
}
