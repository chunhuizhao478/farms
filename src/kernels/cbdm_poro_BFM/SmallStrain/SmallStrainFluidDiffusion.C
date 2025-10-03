//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SmallStrainFluidDiffusion.h"

registerMooseObject("farmsApp", SmallStrainFluidDiffusion);

InputParameters
SmallStrainFluidDiffusion::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Diffusion kernel for pore pressure in poromechanical problems with reference configuration permeability");
  params.addParam<std::string>("base_name", "", "Optional parameter for multiple material systems");
  return params;
}

SmallStrainFluidDiffusion::SmallStrainFluidDiffusion(const InputParameters & parameters)
  : Kernel(parameters),
    _perm_s(getMaterialProperty<Real>("permeability_solid"))
{
}

Real
SmallStrainFluidDiffusion::computeQpResidual()
{

  RankTwoTensor K_s = _perm_s[_qp] * RankTwoTensor::Identity();

  RealVectorValue flux_s = K_s * _grad_u[_qp];

  return _grad_test[_i][_qp] * flux_s;
}

Real
SmallStrainFluidDiffusion::computeQpJacobian()
{
  RankTwoTensor K_s = _perm_s[_qp] * RankTwoTensor::Identity();

  RealVectorValue jac_s = K_s * _grad_phi[_j][_qp];

  return _grad_test[_i][_qp] * jac_s;
}