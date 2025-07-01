//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DamageEvolutionDiffusion.h"

registerMooseObject("farmsApp", DamageEvolutionDiffusion);

InputParameters
DamageEvolutionDiffusion::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("The Laplacian operator ($-nabla cdot nabla u$), with the weak "
                             "form of $(1-B) D nabla alpha cdot nabla phi_i$");
  params.addRequiredCoupledVar("coupled", "The breakage variable");
  return params;
}

DamageEvolutionDiffusion::DamageEvolutionDiffusion(const InputParameters & parameters) 
: Kernel(parameters),
_B_var(coupled("coupled")),
_B(coupledValue("coupled")),
_D_diffusion(getMaterialProperty<Real>("D_diffusion"))
{   
}

Real
DamageEvolutionDiffusion::computeQpResidual()
{
  return (1 - _B[_qp]) * _D_diffusion[_qp] * _grad_u[_qp] * _grad_test[_i][_qp];
}

Real
DamageEvolutionDiffusion::computeQpJacobian()
{
  return (1 - _B[_qp]) * _D_diffusion[_qp] * _grad_phi[_j][_qp] * _grad_test[_i][_qp];
}

Real
DamageEvolutionDiffusion::computeQpOffDiagJacobian(unsigned int jvar)
{
  // d/dB of (1 - B) * D * grad(alpha) · grad(test) = -D * grad(alpha) · grad(test)
  if (jvar == _B_var)
    return -_D_diffusion[_qp] * _grad_u[_qp] * _grad_test[_i][_qp];
  else
    return 0.0;
}