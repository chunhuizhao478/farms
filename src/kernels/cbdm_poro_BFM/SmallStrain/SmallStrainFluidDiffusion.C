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
  params.addClassDescription("Diffusion kernel for pore pressure in poromechanical problems "
                             "with reference configuration permeability and gravity term");
  params.addParam<std::string>("base_name", "", "Optional parameter for multiple material systems");
  params.addParam<RealVectorValue>("gravity_vector", RealVectorValue(0, 0, 0), 
                                   "Gravity vector (default: no gravity)");
  params.addParam<Real>("fluid_density", 0 , "Fluid density");
  return params;
}

SmallStrainFluidDiffusion::SmallStrainFluidDiffusion(const InputParameters & parameters)
  : Kernel(parameters),
    _perm_s(getMaterialProperty<Real>("permeability_solid")),
    _gravity(getParam<RealVectorValue>("gravity_vector")),
    _fluid_density(getParam<Real>("fluid_density"))
{
}

Real
SmallStrainFluidDiffusion::computeQpResidual()
{
  // Permeability tensor (isotropic)
  RankTwoTensor K_s = _perm_s[_qp] * RankTwoTensor::Identity();
  
  // Darcy flux: q = -K * (grad(p) - rho_f * g)
  // Residual: (test, grad(p)) - (test, rho_f * g)
  RealVectorValue flux_s = K_s * (_grad_u[_qp] - _fluid_density * _gravity);
  
  return _grad_test[_i][_qp] * flux_s;
}

Real
SmallStrainFluidDiffusion::computeQpJacobian()
{
  // Permeability tensor (isotropic)
  RankTwoTensor K_s = _perm_s[_qp] * RankTwoTensor::Identity();
  
  // Jacobian only for pressure gradient term (gravity term doesn't depend on u)
  RealVectorValue jac_s = K_s * _grad_phi[_j][_qp];
  
  return _grad_test[_i][_qp] * jac_s;
}