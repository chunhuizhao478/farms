//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "BreakageEvolutionDiffusion.h"

registerMooseObject("farmsApp", BreakageEvolutionDiffusion);

InputParameters
BreakageEvolutionDiffusion::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription(
      "Laplacian contribution for the breakage field: D * grad(B) · grad(phi_i) with the "
      "diffusion coefficient supplied as a material property.");
  return params;
}

BreakageEvolutionDiffusion::BreakageEvolutionDiffusion(const InputParameters & parameters)
  : Kernel(parameters),
    _D_diffusion(getMaterialProperty<Real>("D_diffusion"))
{
}

Real
BreakageEvolutionDiffusion::computeQpResidual()
{
  return _D_diffusion[_qp] * (1.0 - _u[_qp]) * _grad_u[_qp] * _grad_test[_i][_qp];
}

Real
BreakageEvolutionDiffusion::computeQpJacobian()
{
  // d/dB: (1 - B) * grad(B) -> -(grad(B)) + (1 - B) * grad(phi_j)
  const RealVectorValue & grad_phi_j = _grad_phi[_j][_qp];
  const RealVectorValue & grad_u_qp = _grad_u[_qp];

  return _D_diffusion[_qp] *
         ((1.0 - _u[_qp]) * grad_phi_j - _phi[_j][_qp] * grad_u_qp) * _grad_test[_i][_qp];
}
