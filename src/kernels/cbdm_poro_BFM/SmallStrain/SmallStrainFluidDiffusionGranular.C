//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SmallStrainFluidDiffusionGranular.h"

registerMooseObject("farmsApp", SmallStrainFluidDiffusionGranular);

InputParameters
SmallStrainFluidDiffusionGranular::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Diffusion kernel for pore pressure in poromechanical problems with reference configuration permeability");
  params.addParam<std::string>("base_name", "", "Optional parameter for multiple material systems");
  return params;
}

SmallStrainFluidDiffusionGranular::SmallStrainFluidDiffusionGranular(const InputParameters & parameters)
  : Kernel(parameters),
    _perm_g(getMaterialProperty<Real>("permeability_granular"))
{
}

Real
SmallStrainFluidDiffusionGranular::computeQpResidual()
{

  RankTwoTensor K_g = (_perm_g[_qp]) * RankTwoTensor::Identity();
  
  RealVectorValue flux_g = K_g * _grad_u[_qp];

  return _grad_test[_i][_qp] * flux_g;
}

Real
SmallStrainFluidDiffusionGranular::computeQpJacobian()
{
  RankTwoTensor K_g = (_perm_g[_qp]) * RankTwoTensor::Identity();
  
  RealVectorValue jac_flux = K_g * _grad_phi[_j][_qp];

  return _grad_test[_i][_qp] * jac_flux;
}