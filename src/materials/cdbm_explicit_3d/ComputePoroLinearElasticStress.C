//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ComputePoroLinearElasticStress.h"

registerMooseObject("SolidMechanicsApp", ComputePoroLinearElasticStress);

InputParameters
ComputePoroLinearElasticStress::validParams()
{
  InputParameters params = ComputeStressBase::validParams();
  params.addClassDescription("Compute stress using elasticity for small strains");
  params.addRequiredCoupledVar(
    "porepressure", 
    "The pore pressure appropriate for the simulation geometry and coordinate system");
  return params;
}

ComputePoroLinearElasticStress::ComputePoroLinearElasticStress(const InputParameters & parameters)
  : ComputeStressBase(parameters),
    _elasticity_tensor_name(_base_name + "elasticity_tensor"),
    _elasticity_tensor(getMaterialPropertyByName<RankFourTensor>(_elasticity_tensor_name)),
    _pore_pressure(coupledValue("porepressure")),
    _I1(declareProperty<Real>(_base_name + "first_elastic_strain_invariant")),
    _stress_off_diag_jacobian(declareProperty<RankTwoTensor>(_base_name + "stress_off_diag_jacobian")),
    _biot_coeff_eff(getMaterialProperty<Real>(_base_name + "biot_coefficient_effective"))
{
}

void
ComputePoroLinearElasticStress::initialSetup()
{
  // _base_name + "unstabilized_deformation_gradient" is only declared if we're
  // using the Lagrangian kernels.  It's okay to invoke this small strain
  // material if you are using that kernel system and the
  // ComputeLagrangianWrappedStress wrapper
  if (hasBlockMaterialProperty<RankTwoTensor>(_base_name + "strain_increment") &&
      !hasBlockMaterialProperty<RankTwoTensor>(_base_name + "unstabilized_deformation_gradient"))
    mooseError("This linear elastic stress calculation only works for small strains; use "
               "ComputeFiniteStrainElasticStress for simulations using incremental and finite "
               "strains.");
}

void
ComputePoroLinearElasticStress::computeQpStress()
{
  // stress = C * e
  _stress[_qp] = _elasticity_tensor[_qp] * _mechanical_strain[_qp] - _biot_coeff_eff[_qp] * _pore_pressure[_qp] * RankTwoTensor::Identity();

  // Assign value for elastic strain, which is equal to the mechanical strain
  _elastic_strain[_qp] = _mechanical_strain[_qp];

  // Assign value for elastic strain, which is equal to the mechanical strain
  _I1[_qp] = _elastic_strain[_qp].trace();

  // Compute dstress_dstrain
  _Jacobian_mult[_qp] = _elasticity_tensor[_qp];

  // Compute dstress_dpore_pressure
  _stress_off_diag_jacobian[_qp] = -_biot_coeff_eff[_qp] * RankTwoTensor::Identity();
}